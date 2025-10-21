# pySERtool.py - Sequence Pattern Matcher
import os
import argparse
import time
from datetime import datetime
from typing import List, Dict, Any
from openpyxl import Workbook
from Bio import SeqIO
from tqdm import tqdm
from collections import defaultdict
import logging
import sys
import re


# -----------------------------
# Color & Output Configuration
# -----------------------------
COLORS = {
    'success': '\033[38;5;28m',      # Green
    'warn': '\033[38;5;220m',        # Yellow
    'dim': '\033[38;5;245m',         # Gray
    'error': '\033[38;5;196m',       # Red
    'orange': '\033[38;5;208m',      # Orange
    'white': '\033[38;5;231m',       # White
    'end': '\033[0m'
}

def cprint(text: str, color: str = ''):
    """Print colored text to terminal."""
    color_code = COLORS.get(color, '')
    end_code = COLORS['end']
    print(f"{color_code}{text}{end_code}")


# -----------------------------
# Logger Setup
# -----------------------------
def setup_logger(task_name: str, log_dir: str = "./logs"):
    """Setup logger with file and colored console output."""
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"{log_dir}/log_{task_name}_{timestamp}.txt"

    logger = logging.getLogger(task_name)
    if logger.handlers:
        logger.handlers.clear()
    logger.setLevel(logging.INFO)

    # File handler
    fh = logging.FileHandler(log_filename, encoding='utf-8')
    fh.setLevel(logging.INFO)
    file_formatter = logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S')
    fh.setFormatter(file_formatter)
    logger.addHandler(fh)

    # Console handler with colors
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


# -----------------------------
# Core Logic Functions
# -----------------------------
def load_sequences_with_logger(filepath: str, logger: logging.Logger) -> Dict[str, str]:
    """Load sequences from FASTA file using the provided logger."""
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


def calculate_k_b(hit1: int, window1: int, hit2: int = None, window2: int = None) -> (float, float):
    """Calculate k and b for linear formula: window = k * hit + b"""
    if hit1 <= 0:
        raise ValueError("Hit count must be > 0")
    if hit2 is not None and window2 is not None:
        if hit2 <= 0 or hit2 == hit1:
            raise ValueError("Invalid second point")
        k = (window2 - window1) / (hit2 - hit1)
        b = window1 - k * hit1
    else:
        k = window1 / hit1
        b = 0.0
    return k, b


def window_to_hit(window_size: int, k: float, b: float = 0.0) -> int:
    """
    Reverse calculate hit count from window size using linear model: window = k * hit + b
    Solve: hit = (window - b) / k
    """
    hit = (window_size - b) / k
    return max(1, int(hit))


def matches_pattern(pattern: str, text: str) -> bool:
    """Match pattern with wildcard '.'"""
    if len(pattern) != len(text):
        return False
    return all(p == '.' or p == t for p, t in zip(pattern, text))


def scan_sequences_batch(
    sequences: Dict[str, str],
    words: List[str],
    window_size: int,
    hit_count: int,
    slide_step: int = 1
) -> List[Dict]:
    """Scan all sequences for given window/hit condition with custom slide step."""
    hits = []
    for name, seq in sequences.items():
        #⚠️ Warning: slide_step > 1 may miss short motifs
        for i in range(0, len(seq) - window_size + 1, slide_step):
            window = seq[i:i + window_size]
            total_matches = 0
            for word in words:
                word = word.strip()
                if not word:
                    continue
                if '.' in word:
                    for j in range(len(window) - len(word) + 1):
                        if matches_pattern(word, window[j:j + len(word)]):
                            total_matches += 1
                else:
                    total_matches += window.count(word)
            if total_matches >= hit_count:
                hits.append({
                    'gene_name': name,
                    'window_start': i,
                    'window_end': i + window_size,
                    'window_sequence': window,
                    'window_size': window_size,
                    'hit_count': total_matches
                })
    return hits
def is_start_match(seq: str, patterns: List[str]) -> bool:
    escaped = [re.escape(p).replace('\\.', '.') for p in patterns]
    return any(re.match(f"^{p}", seq) for p in escaped)


def is_end_match(seq: str, patterns: List[str]) -> bool:
    escaped = [re.escape(p).replace('\\.', '.') for p in patterns]
    return any(re.match(f".*{p}$", seq) for p in escaped)

def save_to_excel(unique_sequences: Dict[str, Any], excel_file: str, pattern_length: int, all_words: List[str]):
    wb = Workbook()
    if wb.active:
        wb.remove(wb.active)
    ws = wb.create_sheet("Unique Sequences")

    headers = [
        'Item titles',
        'Sequence',
        'Matched Windows',
        'Configurations',
        'Max Diff Info',
        'Score'
    ]
    ws.append(headers)

    num_dots = all_words[0].count('.')

    all_rows = []

    for gene, data in unique_sequences.items():
        matched_windows = data['Matched Windows']
        if not matched_windows:
            continue

        window_scores = []
        for w in matched_windows:
            raw_net_diff = 2 * pattern_length * w['Hit Count'] - w['Window Size']
            penalized_net_diff = raw_net_diff - num_dots * w['Window Size']
            window_scores.append((penalized_net_diff, w))

        max_net_diff = max(window_scores, key=lambda x: x[0])[0]
        best_raw_windows = [w for score, w in window_scores if score == max_net_diff]

        trimmed_candidates = []

        for win in best_raw_windows:
            seq = win['Window Sequence']
            start = win['Window Start']
            end = win['Window End'] - 1
            current_seq = seq
            current_start = start
            current_end = end
            trim_count = 0

            while len(current_seq) > 0:
                if is_start_match(current_seq, all_words):
                    break
                current_seq = current_seq[1:]
                current_start += 1
                trim_count += 1

            while len(current_seq) > 0:
                if is_end_match(current_seq, all_words):
                    break
                current_seq = current_seq[:-1]
                current_end -= 1
                trim_count += 1

            if len(current_seq) == 0:
                continue

            adjusted_score = max_net_diff + trim_count
            info_str = f"[{current_start}-{current_end}] {adjusted_score:+d}: {current_seq}"
            trimmed_candidates.append({
                'score': adjusted_score,
                'info': info_str,
                'start': current_start,
                'end': current_end,
                'sequence': current_seq
            })

        if not trimmed_candidates:
            max_diff_info = "N/A"
            final_score = 0
        else:
            max_score = max(c['score'] for c in trimmed_candidates)
            top_candidates = [c for c in trimmed_candidates if c['score'] == max_score]
            seen = set()
            unique_top = []
            for c in top_candidates:
                key = (c['start'], c['end'], c['sequence'])
                if key not in seen:
                    seen.add(key)
                    unique_top.append(c)
            unique_top.sort(key=lambda x: x['start'])
            max_diff_info = "; ".join(c['info'] for c in unique_top)
            final_score = max_score

        window_str = "; ".join(
            f"{w['Window Start']}-{w['Window End']-1}:{w['Window Sequence']}"
            for w in sorted(matched_windows, key=lambda x: x['Window Start'])
        )

        configs = {(w['Window Size'], w['Hit Count']) for w in matched_windows}
        config_str = ", ".join(f"(W{w},H{h})" for w, h in sorted(configs))

        all_rows.append([
            gene,
            data['Sequence'],
            window_str,
            config_str,
            max_diff_info,
            final_score
        ])

    all_rows.sort(key=lambda x: (-x[5], x[0]))

    for row in all_rows:
        ws.append(row)

    wb.save(excel_file)
    cprint(f"✅ Excel saved: {excel_file}", 'success')

# -----------------------------
# Task Processing
# -----------------------------
def process_task(
    sequences: Dict[str, str],
    task_words: str,
    start_window: int,
    end_window: int,
    window_step: int,
    slide_step: int,
    mode: str,
    points: List[int],
    output_dir: str,
    input_filename: str
) -> Dict[str, Any]:
    """Process a single task with progress bar and logging."""
    words = [w.strip() for w in task_words.split('-') if w.strip()]
    if not words:
        raise ValueError("No valid words in task")
    pattern_length = len(words[0])
    input_base = os.path.splitext(os.path.basename(input_filename))[0]
    excel_file = os.path.join(output_dir, f"{input_base}_{task_words}.xlsx")

    logger, log_file = setup_logger(task_words)
    points_str = ",".join(map(str, points))
    cprint(f"Task: '{task_words}' | Window: {start_window}-{end_window} (wstep={window_step}) | Slide Step: {slide_step} | Mode: {mode} | Points: {points_str}", 'orange')
    if slide_step > 1:
        logger.warning(f"Using slide_step={slide_step} may miss short motifs. Use 1 for full accuracy.")
    configs = []
    if mode == 'formula':
        h1, w1 = points[0], points[1]
        if len(points) == 4:
            h2, w2 = points[2], points[3]
            k, b = calculate_k_b(h1, w1, h2, w2)
        else:
            k, b = calculate_k_b(h1, w1)
        configs = [
            (w, window_to_hit(w, k, b))
            for w in range(start_window, end_window + 1, window_step)]
        logger.info(f"Formula mode: k={k:.6f}, b={b:.6f}")
    else:
        if len(points) % 2 != 0:
            raise ValueError("In 'fix' mode, --point must have even number of integers")
        configs = [(points[i+1], points[i]) for i in range(0, len(points), 2)]
        configs = [(w, h) for w, h in configs if start_window <= w <= end_window]

    logger.info(f"Generated {len(configs)} raw configurations")

    # Max-Window per Hit Strategy (HMS)
    max_window_per_hit = defaultdict(int)
    for w, h in configs:
        if w > max_window_per_hit[h]:
            max_window_per_hit[h] = w

    final_configs = [(max_window_per_hit[h], h) for h in sorted(max_window_per_hit.keys())]
    final_configs.sort(key=lambda x: x[0])

    # Log HMS strategy
    items = []
    for h in sorted(max_window_per_hit.keys()):
        max_w = max_window_per_hit[h]
        if mode == 'formula':
            covered_ws = [w for w in range(start_window, end_window + 1, window_step) if window_to_hit(w, k, b) == h]
        else:
            covered_ws = [w for w, hit in configs if hit == h]
        covered_str = ",".join([f"W{w}" for w in sorted(covered_ws)])
        items.append(f"W{max_w}:H{h:<2} ({covered_str})")

    logger.info("Max-Window Coverage Strategy (HMS):")
    for i in range(0, len(items), 4):
        line = "  ".join(items[i:i+4])
        logger.info(f"   {line}")

    logger.info(f"Scanning {len(final_configs)} optimized configs (was {len(configs)})")

    # Run scanning
    unique_sequences = {}
    task_start = time.time()
    start_time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    logger.info(f"[{start_time_str}] Starting task: {task_words}")

    for w, h in tqdm(
        final_configs,
        desc=f"🔍 Scanning '{task_words}'",
        unit="config",
        ncols=100,
        colour='green',
        leave=True,
        file=sys.stdout
    ):
        hits = scan_sequences_batch(sequences, words, w, h, slide_step)
        for hit in hits:
            gene_name = hit['gene_name']
            if gene_name not in unique_sequences:
                unique_sequences[gene_name] = {
                    'Sequence': sequences[gene_name],
                    'Matched Windows': []}
            unique_sequences[hit['gene_name']]['Matched Windows'].append({
                'Window Start': hit['window_start']+1,
                'Window End': hit['window_end']+1,
                'Window Sequence': hit['window_sequence'],
                'Window Size': hit['window_size'],
                'Hit Count': hit['hit_count']
            })

    task_time = time.time() - task_start
    end_time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    genes_found = len(unique_sequences)

    save_to_excel(unique_sequences, excel_file, pattern_length=pattern_length, all_words=words)

    logger.info(f"[{end_time_str}] Task '{task_words}' completed!")
    cprint(f"✅ Found {genes_found} unique genes | 💾 Saved: {os.path.basename(excel_file)}", 'white')
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


# -----------------------------
# Main Entry Point
# -----------------------------
def main():
    total_start = time.time()
    parser = argparse.ArgumentParser(description="Gene Pattern Matcher with Wildcard and Step Support")
    parser.add_argument('tasks', type=str, help="Tasks: e.g., AT-CG,GCTA or D.E,F.G")
    parser.add_argument('--input', type=str, default='./Database/normal.fasta', help="FASTA input file path")
    parser.add_argument('--start', type=int, default=20, help="Start window size")
    parser.add_argument('--end', type=int, default=100, help="End window size")
    parser.add_argument('--step', type=int, default=1, help="Stride step for sliding window (how many bases to move right)")
    parser.add_argument('--wstep', type=int, default=1, help="Step size for window size increment (default: 1)")
    parser.add_argument('--output', type=str, default='./Result', help="Output directory")
    parser.add_argument('--point', nargs='+', type=int, help="Calibration points: H1 W1 [H2 W2] for formula, or list of H-W pairs for fix")
    parser.add_argument('--mode', choices=['formula', 'fix'], default='formula', help="Scanning mode")

    args = parser.parse_args()

    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")
    os.makedirs(args.output, exist_ok=True)

    cprint(f"\n🚀 Starting {len(args.tasks.split(','))} tasks...", 'warn')

    first_task_name = args.tasks.split(',')[0].strip() or "GeneMatcher"
    logger, _ = setup_logger(first_task_name)
    sequences = load_sequences_with_logger(args.input, logger)

    if args.mode == 'formula':
        if args.point is None:
            points = [18, 20, 30, 50]
            logger.warning(f"Using default calibration points: {points}")
            cprint(f"⚠️  Using default points: {points}", 'warn')
        elif len(args.point) == 2 or len(args.point) == 4:
            points = args.point
        else:
            raise ValueError("In 'formula' mode, --point must have 2 or 4 integers")
    else:
        if args.point is None or len(args.point) % 2 != 0:
            raise ValueError("--point must provide an even number of integers in 'fix' mode")
        points = args.point

    results = []
    for task in args.tasks.split(','):
        task = task.strip()
        if not task:
            continue
        result = process_task(
            sequences=sequences,
            task_words=task,
            start_window=args.start,
            end_window=args.end,
            window_step=args.wstep,
            slide_step=args.step,
            mode=args.mode,
            points=points,
            output_dir=args.output,
            input_filename=args.input
        )
        results.append(result)

    total_duration = time.time() - total_start

    cprint("=" * 60, 'dim')
    cprint("🎉 ALL TASKS COMPLETED", 'success')
    cprint(f"⏱  Total Time: {total_duration:.2f}s", 'success')
    cprint(f"📁 Output Directory: {args.output}", 'dim')
    cprint(f"📏 Window Range: {args.start}-{args.end} (wstep={args.wstep})", 'dim')
    cprint(f"🧩 Slide Step: {args.step}", 'dim')
    cprint(f"⚙️  Mode: {args.mode.upper()}", 'dim')
    cprint("=" * 60, 'dim')

    cprint("📊 Summary:", 'warn')
    for r in results:
        cprint(f"  • {r['task']:12} → {r['genes_found']:3} genes | {r['time']:6.2f}s", 'success')


if __name__ == "__main__":
    main()
