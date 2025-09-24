# SERtool.py
import os
import argparse
import time
from datetime import datetime
from typing import List, Dict, Any
from openpyxl import Workbook
from Bio import SeqIO
import ser_tool
from tqdm import tqdm
from collections import defaultdict
import logging
import sys


COLORS = {
    'success': '\033[38;5;28m',
    'warn': '\033[38;5;220m',
    'dim': '\033[38;5;245m',
    'error': '\033[38;5;196m',
    'orange': '\033[38;5;208m',
    'white': '\033[38;5;231m',
    'end': '\033[0m'
}

def cprint(text: str, color: str = ''):
    color_code = COLORS.get(color, '')
    end_code = COLORS['end']
    print(f"{color_code}{text}{end_code}")


def setup_logger(task_name: str, log_dir: str):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"{log_dir}/log_{task_name}_{timestamp}.txt"

    logger = logging.getLogger(task_name + "_" + timestamp)
    if logger.handlers:
        logger.handlers.clear()
    logger.setLevel(logging.INFO)

    fh = logging.FileHandler(log_filename, encoding='utf-8')
    fh.setLevel(logging.INFO)
    file_formatter = logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S')
    fh.setFormatter(file_formatter)
    logger.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    class DarkFormatter(logging.Formatter):
        def format(self, record):
            level_colors = {
                logging.INFO: COLORS['success'],
                logging.WARNING: COLORS['warn'],
                logging.ERROR: COLORS['error'],
            }
            end = COLORS['end']
            color = level_colors.get(record.levelno, '')
            if not record.msg.startswith(('\033[', end)):
                record.msg = f"{color}{record.msg}{end}"
            return super().format(record)

    ch.setFormatter(DarkFormatter('%(message)s', datefmt='%H:%M:%S'))
    logger.addHandler(ch)

    return logger, log_filename


def load_sequences(filepath: str, logger: logging.Logger) -> Dict[str, str]:
    logger.info(f"Loading sequences from {filepath}...")
    start = time.time()
    try:
        sequences = {record.description: str(record.seq) for record in SeqIO.parse(filepath, "fasta")}
    except Exception as e:
        logger.error(f"Failed to read FASTA file: {e}")
        raise
    duration = time.time() - start
    logger.info(f"Loaded {len(sequences)} sequences in {duration:.2f}s")
    return sequences


def calculate_hit(window_size: int, k: float, b: float) -> int:
    hit = int((window_size - b) / k)
    return max(1, hit)


def save_to_excel(unique_sequences: Dict[str, Any], excel_file: str, pattern_length: int):
    wb = Workbook()
    wb.remove(wb.active)
    ws = wb.create_sheet("Unique Sequences")
    headers = ['Gene Name', 'Sequence', 'Matched Windows', 'Configurations', 'Max Diff Window (Match - Unmatch)']
    ws.append(headers)

    for gene, data in unique_sequences.items():
        matched = data['Matched Windows']
        window_str = "; ".join(f"{w['Window Start']}-{w['Window End']}:{w['Window Sequence']}" for w in matched)
        configs = {(w['Window Size'], w['Hit Count']) for w in matched}
        config_str = ", ".join(f"(W{w},H{h})" for w, h in sorted(configs))
        best = max(matched, key=lambda x: 2 * pattern_length * x['Hit Count'] - x['Window Size'])
        matched_chars = pattern_length * best['Hit Count']
        unmatched_chars = best['Window Size'] - matched_chars
        net_diff = matched_chars - unmatched_chars
        max_diff_info = (
            f"[{best['Window Start']}-{best['Window End']}] "
            f"{net_diff:+d} ({matched_chars}M-{unmatched_chars}U): {best['Window Sequence']}")
        ws.append([gene, data['Sequence'], window_str, config_str, max_diff_info])
    wb.save(excel_file)


def process_task(
    sequences: Dict[str, str],
    task_words: str,
    start_window: int,
    end_window: int,
    step: int,
    mode: str,
    points: List[int],
    shared_output_dir: str,
    input_filename: str
) -> Dict[str, Any]:
    words = [w.strip() for w in task_words.split('-') if w.strip()]
    
    if not words:
        raise ValueError("No valid task words provided")

    pattern_length = len(words[0])
    lengths = [len(w) for w in words]
    if not all(l == pattern_length for l in lengths):
        cprint(f"Warning: Task words have different lengths: {dict(zip(words, lengths))}. Using {pattern_length}.", 'warn')
    logger, log_file = setup_logger(task_words, log_dir=shared_output_dir)

    task_safe_name = task_words.replace("-", "_").replace(" ", "")
    excel_file = os.path.join(shared_output_dir, f"{task_safe_name}.xlsx")

    points_str = ",".join(map(str, points))
    cprint(f"Task: '{task_words}' | Window: {start_window}-{end_window} (step={step}) | Points: {points_str}", 'orange')

    logger.info(f"   **Configuration:")
    logger.info(f"   Mode       : {mode}")
    logger.info(f"   Window     : {start_window} - {end_window} (step={step})")
    logger.info(f"   Points     : [{points_str}]")
    logger.info(f"   Task Words : {task_words}")
    logger.info(f"   Input File : {input_filename}")
    logger.info(f"   Seq Count  : {len(sequences)}")

    configs = []
    if mode == 'formula':
        h1, w1 = points[0], points[1]
        if len(points) == 4:
            h2, w2 = points[2], points[3]
            k, b = ser_tool.calculate_k_b(h1, w1, h2, w2)
        else:
            k, b = ser_tool.calculate_k_b(h1, w1, None, None)
        configs = [
            (w, calculate_hit(w, k, b))
            for w in range(start_window, end_window + 1, step)
        ]
    else:
        configs = [(points[i+1], points[i]) for i in range(0, len(points), 2)]

    logger.info(f"Original {len(configs)} configurations")

    # Max-Window per Hit Strategy
    max_window_per_hit = defaultdict(int)
    for w, h in configs:
        if w > max_window_per_hit[h]:
            max_window_per_hit[h] = w

    final_configs = [(max_window_per_hit[h], h) for h in sorted(max_window_per_hit.keys())]
    final_configs.sort(key=lambda x: x[0])

    items = []
    for h in sorted(max_window_per_hit.keys()):
        max_w = max_window_per_hit[h]
        if mode == 'formula':
            covered_ws = [w for w in range(start_window, end_window + 1, step)
                          if calculate_hit(w, k, b) == h]
        else:
            covered_ws = [w for w, hit in configs if hit == h]
        covered_str = ",".join([f"W{w}" for w in sorted(covered_ws)])
        items.append(f"W{max_w}:H{h:<2} ({covered_str})")

    logger.info("Max-Window Coverage Strategy (HMS):")
    for i in range(0, len(items), 4):
        line = "  ".join(items[i:i+4])
        logger.info(f"   {line}")

    logger.info(f"Scanning {len(final_configs)} configs (was {len(configs)})")

    unique_sequences = {}
    task_start = time.time()
    start_time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    logger.info(f"[{start_time_str}] Starting task: {task_words}")

    for w, h in tqdm(
        final_configs,
        desc=f"Scanning '{task_words}'",
        unit="config",
        ncols=100,
        colour='green',
        leave=True,
        file=sys.stdout
    ):
        hits = ser_tool.scan_sequences_batch(sequences, words, w, h)
        for hit in hits:
            if hit.gene_name not in unique_sequences:
                unique_sequences[hit.gene_name] = {
                    'Sequence': sequences[hit.gene_name],
                    'Matched Windows': []
                }
            unique_sequences[hit.gene_name]['Matched Windows'].append({
                'Window Start': hit.window_start,
                'Window End': hit.window_end,
                'Window Sequence': hit.window_sequence,
                'Window Size': hit.window_size,
                'Hit Count': hit.hit_count
            })

    task_time = time.time() - task_start
    end_time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    genes_found = len(unique_sequences)

    save_to_excel(unique_sequences, excel_file, pattern_length=pattern_length)

    logger.info(f"[{end_time_str}] Task '{task_words}' completed!")
    cprint(f"Found {genes_found} unique genes | 💾 Saved: {os.path.basename(excel_file)}", 'white')
    logger.info(f"Duration: {task_time:.2f}s")
    logger.info(f"Full log saved to: {log_file}")

    return {
        'task': task_words,
        'time': task_time,
        'genes_found': genes_found,
        'excel': excel_file,
        'start_time': start_time_str,
        'end_time': end_time_str,
        'log_file': log_file
    }


def main():
    total_start = time.time()
    parser = argparse.ArgumentParser(description="Sequence scanner with compact HMS output")
    parser.add_argument('tasks', type=str, help="Task words, e.g., AB-CD,TG-CA")
    parser.add_argument('--input', type=str, default='./Database/normal.fasta', help="Input FASTA file")
    parser.add_argument('--start', type=int, default=20, help="Start window size")
    parser.add_argument('--end', type=int, default=100, help="End window size")
    parser.add_argument('--step', type=int, default=1, help="Step size for window")
    parser.add_argument('--output', type=str, default='./Result', help="Output directory")
    parser.add_argument('--point', nargs='+', type=int, help="Calibration points")
    parser.add_argument('--mode', choices=['formula', 'fix'], default='formula', help="Scanning mode")

    args = parser.parse_args()

    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")
    os.makedirs(args.output, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    point_str = "-".join(map(str, args.point)) if args.point else "default"
    folder_name = f"{timestamp}_W{args.start}-{args.end}_S{args.step}_{args.mode}_P{point_str}"
    shared_output_dir = os.path.join(args.output, folder_name)
    os.makedirs(shared_output_dir, exist_ok=True)

    cprint(f"\nOutput Directory: {shared_output_dir}", 'dim')
    cprint(f"Tasks: {args.tasks}", 'warn')

    temp_logger, _ = setup_logger("LOAD", log_dir=shared_output_dir)
    sequences = load_sequences(args.input, temp_logger)

    if args.mode == 'formula':
        if args.point is None:
            points = [18, 20, 30, 50]
            cprint(f"Using default points: {points}", 'warn')
        else:
            points = args.point
    else:
        if args.point is None or len(args.point) % 2 != 0:
            raise ValueError("--point must be even number of ints in 'fix' mode")
        points = args.point

    results = []
    for task in args.tasks.split(','):
        task = task.strip()
        if not task:
            continue
        result = process_task(
            sequences, task, args.start, args.end, args.step,
            args.mode, points, shared_output_dir, args.input
        )
        results.append(result)

    total_duration = time.time() - total_start

    cprint("=" * 60, 'dim')
    cprint("ALL TASKS COMPLETED", 'success')
    cprint(f"Total Time: {total_duration:.2f}s", 'success')
    cprint(f"Output Folder: {folder_name}", 'dim')
    cprint(f"Root Path: {shared_output_dir}", 'dim')
    cprint("=" * 60, 'dim')

    cprint("Summary:", 'warn')
    for r in results:
        cprint(f"  • {r['task']:12} → {r['genes_found']:3} genes | {r['time']:6.2f}s", 'success')

if __name__ == "__main__":

    main()
