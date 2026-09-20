"""
Build marker-specific BLAST databases from DNABarcode records.

This command groups DNA barcodes by marker family and creates separate
BLAST databases for each marker type (COI, ITS, 18S, matK, rbcL, etc.).

Usage:
    python manage.py build_marker_blast_db
"""
import os
import subprocess
import tempfile
from pathlib import Path

from django.core.management.base import BaseCommand
import pests_project.settings as app_settings

from quarantine_pests.models import DNABarcode

# 定义Marker家族成员映射
MARKER_FAMILY_MEMBERS = {
    'COI': ['COI-5P', 'COI-3P', 'COII', 'COXIII', 'COI-PSEUDO', 'COI-5PNMT1', 'COII-COI'],
    'ITS': ['ITS', 'ITS1', 'ITS2', '5-8S'],
    '18S': ['18S', '18S-5P', '18S-3P'],
    'matK': ['matK', 'matK-like'],
    'rbcL': ['rbcL', 'rbcLa', 'rbcL-like'],
    '16S': ['16S'],
    'CYTB': ['CYTB'],
    '28S': ['28S', '28S-D2', '28S-D2-D3', '28S-D1-D2', '28S-D3'],
    'ND': ['ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5-0', 'ND6', 'MT-NADH5'],
    'trn': ['trnH-psbA', 'trnL', 'trnL-F', 'trnD-trnY-trnE', 'cp-trnL-intron', 'atpF-atpH',
            'atpF-intron', 'atpB-rbcL', 'ndhK-ndhC', 'psbK-psbI', 'psbE-petL',
            'rpL32-trnL', 'rps16-trnK'],
    'OTHER': None,  # 动态计算：排除所有已知家族
}


class Command(BaseCommand):
    help = 'Build marker-specific BLAST databases from DNABarcode records'

    def handle(self, *args, **options):
        # Define marker family to DB path mapping
        marker_db_config = {
            'COI': app_settings.BLAST_DB_COI,
            'ITS': app_settings.BLAST_DB_ITS,
            '18S': app_settings.BLAST_DB_18S,
            'matK': app_settings.BLAST_DB_matK,
            'rbcL': app_settings.BLAST_DB_rbcL,
            '16S': app_settings.BLAST_DB_16S,
            'CYTB': app_settings.BLAST_DB_CYTB,
            '28S': app_settings.BLAST_DB_28S,
            'ND': app_settings.BLAST_DB_ND,
            'trn': app_settings.BLAST_DB_trn,
            'OTHER': app_settings.BLAST_DB_OTHER,
        }

        makeblastdb_path = getattr(app_settings, 'BLAST_MAKEBLASTDB_BIN_PATH', None)
        if not makeblastdb_path:
            self.stderr.write(self.style.ERROR(
                'BLAST_MAKEBLASTDB_BIN_PATH not configured in settings'
            ))
            return

        if not os.path.exists(makeblastdb_path):
            self.stderr.write(self.style.ERROR(
                f'makeblastdb not found at: {makeblastdb_path}'
            ))
            return

        fasta_paths = []  # 收集所有FASTA文件路径，用于构建ALL数据库

        for marker, db_path in marker_db_config.items():
            self.stdout.write(f'\n{"="*60}')
            self.stdout.write(f'Processing marker family: {marker}')
            self.stdout.write(f'{"="*60}')

            # Get barcodes for this marker family
            family_members = MARKER_FAMILY_MEMBERS.get(marker)
            if family_members is None:
                # OTHER: 排除所有已知家族的成员
                all_known_markers = []
                for members in MARKER_FAMILY_MEMBERS.values():
                    if members:
                        all_known_markers.extend(members)
                barcodes = DNABarcode.objects.exclude(
                    marker_code__in=all_known_markers
                ).exclude(sequence__isnull=True).exclude(sequence='')
            else:
                barcodes = DNABarcode.objects.filter(
                    marker_code__in=family_members
                ).exclude(sequence__isnull=True).exclude(sequence='')

            count = barcodes.count()
            self.stdout.write(f'Found {count} barcodes for marker family {marker}')
            if family_members:
                self.stdout.write(f'  Included markers: {", ".join(family_members)}')

            if count == 0:
                self.stdout.write(self.style.WARNING(
                    f'No barcodes found for {marker}, skipping database build'
                ))
                continue

            # Create DB directory if it doesn't exist
            db_path = Path(db_path)
            db_path.mkdir(parents=True, exist_ok=True)

            # Generate FASTA file
            fasta_path = db_path / 'barcodes.fasta'
            self._write_fasta(barcodes, fasta_path, marker)
            self.stdout.write(f'Written FASTA to: {fasta_path}')

            # 收集FASTA文件路径
            fasta_paths.append(fasta_path)

            # Build BLAST database
            self._build_blast_db(makeblastdb_path, fasta_path, db_path, marker)

        # Build ALL database (combined all markers)
        self._build_all_database(makeblastdb_path, fasta_paths)

        self.stdout.write(self.style.SUCCESS('\nAll marker-specific BLAST databases built successfully!'))

    def _write_fasta(self, barcodes, fasta_path, marker):
        """Write barcodes to FASTA file, using process_id as header."""
        import re
        valid_nucleotide_chars = set('ATCGN')  # 只允许有效的核苷酸字符

        # 定义最小长度要求
        # 所有 marker 统一采用 200 bp 的最小长度标准
        min_length = 200

        with open(fasta_path, 'w', encoding='utf-8') as f:
            seen_ids = set()  # 用于跟踪已使用的 ID，处理重复的 process_id
            skipped_invalid = 0
            skipped_duplicate = 0
            skipped_short = 0

            for barcode in barcodes.iterator():
                # Clean sequence - remove whitespace and convert to uppercase
                sequence = ''.join(barcode.sequence.split()).upper()
                if not sequence:
                    continue

                # 验证序列只包含有效的核苷酸字符
                seq_chars = set(sequence)
                if not seq_chars.issubset(valid_nucleotide_chars):
                    # 跳过包含非法字符的序列
                    skipped_invalid += 1
                    continue

                # 验证序列长度
                if len(sequence) < min_length:
                    skipped_short += 1
                    continue

                # 处理重复的 process_id：跳过已处理过的
                if barcode.process_id in seen_ids:
                    skipped_duplicate += 1
                    continue
                seen_ids.add(barcode.process_id)

                # Write FASTA format: >process_id\nsequence\n
                f.write(f'>{barcode.process_id}\n')
                # Write sequence in lines of 80 characters
                for i in range(0, len(sequence), 80):
                    f.write(sequence[i:i+80] + '\n')

            if skipped_invalid > 0 or skipped_duplicate > 0 or skipped_short > 0:
                self.stdout.write(f'  Skipped {skipped_invalid} invalid sequences, {skipped_short} too short (<{min_length}bp), {skipped_duplicate} duplicates')

    def _build_blast_db(self, makeblastdb_path, fasta_path, db_path, marker):
        """Run makeblastdb to create the BLAST database."""
        self.stdout.write(f'Building BLAST database for {marker}...')

        try:
            result = subprocess.run(
                [
                    makeblastdb_path,
                    '-in', str(fasta_path),
                    '-parse_seqids',
                    '-dbtype', 'nucl',
                    '-out', str(db_path / 'barcodes'),
                ],
                capture_output=True,
                text=True,
                timeout=300
            )

            if result.returncode == 0:
                self.stdout.write(self.style.SUCCESS(
                    f'  Successfully created BLAST database for {marker}'
                ))
            else:
                self.stderr.write(self.style.ERROR(
                    f'  Failed to create BLAST database for {marker}: {result.stderr}'
                ))

        except subprocess.TimeoutExpired:
            self.stderr.write(self.style.ERROR(
                f'  Timeout while building database for {marker}'
            ))
        except Exception as e:
            self.stderr.write(self.style.ERROR(
                f'  Error building database for {marker}: {str(e)}'
            ))

    def _build_all_database(self, makeblastdb_path, fasta_paths):
        """Build a combined BLAST database from all marker-specific FASTA files."""
        if not fasta_paths:
            self.stdout.write(self.style.WARNING('No FASTA files to combine for ALL database'))
            return

        self.stdout.write(f'\n{"="*60}')
        self.stdout.write('Building combined ALL database')
        self.stdout.write(f'{"="*60}')

        # Create ALL database directory
        all_db_path = Path(app_settings.BLAST_DB_ALL)
        all_db_path.mkdir(parents=True, exist_ok=True)

        # Combine all FASTA files with cross-database deduplication
        combined_fasta_path = all_db_path / 'barcodes.fasta'
        total_records = 0
        skipped_duplicates = 0
        seen_ids = set()  # 全局去重：跨所有marker数据库

        with open(combined_fasta_path, 'w', encoding='utf-8') as out_f:
            for fasta_path in fasta_paths:
                marker_name = fasta_path.parent.name
                if fasta_path.exists():
                    with open(fasta_path, 'r', encoding='utf-8') as in_f:
                        content = in_f.read()
                        # 按记录处理，保留第一次出现的ID
                        records = content.split('>')
                        for record in records:
                            if not record.strip():
                                continue
                            lines = record.strip().split('\n')
                            record_id = lines[0].split('\s')[0]  # 取ID部分，去除可能的长度信息
                            if record_id in seen_ids:
                                skipped_duplicates += 1
                                continue
                            seen_ids.add(record_id)
                            out_f.write('>' + record + '\n')
                            total_records += 1
                    self.stdout.write(f'  Added {marker_name}')

        if skipped_duplicates > 0:
            self.stdout.write(f'  Skipped {skipped_duplicates} duplicate records across databases')
        self.stdout.write(f'Combined {len(fasta_paths)} FASTA files with {total_records} total unique records')

        # Build BLAST database for ALL
        self._build_blast_db(makeblastdb_path, combined_fasta_path, all_db_path, 'ALL')