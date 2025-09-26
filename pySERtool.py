import os
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from openpyxl import Workbook
from Bio import SeqIO


def load_sequences(filepath):
    """Return {full description line: sequence}"""
    return {record.description: str(record.seq) 
            for record in SeqIO.parse(filepath, "fasta")}


def calculate_k_b(hit1, window1, hit2=None, window2=None):
    """
    Fit the linear formula: window = k * hit + b
    Support single point (through origin) or two points mode.
    """
    if hit1 <= 0:
        raise ValueError("Hit count must be greater than 0")

    if hit2 is not None and window2 is not None:
        if hit2 <= 0:
            raise ValueError("Hit count must be greater than 0")
        if hit2 == hit1:
            raise ValueError("Hit counts of two points cannot be the same (k cannot be calculated)")
        k = (window2 - window1) / (hit2 - hit1)
        b = window1 - k * hit1
    else:
        # Single point mode: assume through origin (0,0)
        k = window1 / hit1
        b = 0.0
    return k, b


def calculate_hit(window, k, b=0.0):
    """
    Reverse calculate hit count from linear formula: window = k * hit + b
    Return integer hit (floor, but at least 1)
    """
    if k <= 0:
        raise ValueError("k must be greater than 0")
    value = (window - b) / k
    if value < 1:
        return 1
    return int(value)


def matches_pattern(pattern, text):
    """Pattern matching with wildcard '.'"""
    if len(pattern) != len(text):
        return False
    return all(p == '.' or p == t for p, t in zip(pattern, text))


def process_sequences_with_wind_hit(sequences, words, wind, hit):
    """
    Search for windows meeting hit condition in all sequences
    Return: all hits + unique sequences (by Gene Name)
    """
    all_hits = []
    unique_sequences = {}

    for name, seq in sequences.items():
        for i in range(len(seq) - wind + 1):
            window = seq[i:i + wind]
            total_matches = 0

            for word in words:
                word = word.strip()
                if not word:
                    continue
                if '.' in word:
                    # wildcard matching
                    for j in range(len(window) - len(word) + 1):
                        if matches_pattern(word, window[j:j + len(word)]):
                            total_matches += 1
                else:
                    total_matches += window.count(word)

            if total_matches >= hit:
                # record hit details
                hit_info = {
                    'Gene Name': name,
                    'Sequence': seq,
                    'Window Size': wind,
                    'Hit Count': hit,
                    'Window Start': i,
                    'Window End': i + wind,
                    'Window Sequence': window
                }
                all_hits.append(hit_info)

                # record for deduplication
                if name not in unique_sequences:
                    unique_sequences[name] = {
                        'Gene Name': name,
                        'Sequence': seq,
                        'Window Size': wind,
                        'Hit Count': hit
                    }

    return all_hits, list(unique_sequences.values())


def handle_task(sequences, task_words, start_window, end_window, mode, points, input_filename, output_dir, outputtype, k_val=None, b_val=None):
    """
    Process a single task and generate independent Excel file
    Modified: in simple mode, last column shows max hit/window ratio with corresponding window sequence
    """
    task_name = '_'.join(task_words.replace('.', '_').split('-'))
    input_base = os.path.splitext(os.path.basename(input_filename))[0]
    excel_file = os.path.join(output_dir, f"{input_base}_{task_words}.xlsx")

    wb = Workbook()
    wb.remove(wb.active)  # remove default sheet

    all_hits_entries = []
    unique_sequences_entries = {}

    if mode == 'formula':
        if k_val is None or b_val is None:
            k, b = calculate_k_b(*points)
        else:
            k, b = k_val, b_val
        print(f"Task '{task_words}': using k={k:.6f}, b={b:.6f}")

        config_list = {}
        for w in range(start_window, end_window + 1):
            try:
                h = calculate_hit(w, k, b)
                if h not in config_list or w > config_list[h]:
                    config_list[h] = w
            except ValueError:
                continue
        config_list = [(config_list[h], h) for h in sorted(config_list.keys())]

    elif mode == 'fix':
        if len(points) % 2 != 0:
            raise ValueError("In fix mode, number of --point arguments must be even (x1 y1 x2 y2 ...)")
        config_list = [(points[i], points[i + 1]) for i in range(0, len(points), 2)]
    else:
        raise ValueError("mode must be 'formula' or 'fix'")

    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(process_sequences_with_wind_hit, sequences, task_words.split('-'), w, h)
            for w, h in config_list
        ]

        for future in as_completed(futures):
            hits, uniques = future.result()
            all_hits_entries.extend(hits)

            for hit_info in hits:
                gene = hit_info['Gene Name']
                if gene not in unique_sequences_entries:
                    unique_sequences_entries[gene] = {
                        'Gene Name': gene,
                        'Sequence': hit_info['Sequence'],
                        'Matched Windows': []
                    }
                window_info = {
                    'Window Start': hit_info['Window Start'],
                    'Window End': hit_info['Window End'],
                    'Window Sequence': hit_info['Window Sequence'],
                    'Window Size': hit_info['Window Size'],
                    'Hit Count': hit_info['Hit Count']
                }
                unique_sequences_entries[gene]['Matched Windows'].append(window_info)

    # Write to Excel
    if outputtype == 'all':
        ws_all = wb.create_sheet("All Hits")
        headers_all = ['Gene Name', 'Sequence', 'Window Size', 'Hit Count', 'Window Start', 'Window End', 'Window Sequence']
        ws_all.append(headers_all)
        for row in all_hits_entries:
            ws_all.append([row[h] for h in headers_all])

    ws_unique = wb.create_sheet("Unique Sequences")
    headers_unique = ['Gene Name', 'Sequence', 'Matched Windows', 'Window Size', 'Max Rate Window']
    ws_unique.append(headers_unique)

    for gene_data in unique_sequences_entries.values():
        matched_windows = gene_data['Matched Windows']

        # Format all matched windows
        window_str = "; ".join(
            f"{w['Window Start']}-{w['Window End']}:{w['Window Sequence']}"
            for w in matched_windows
        )

        # List all (window, hit) configurations
        configs = {(w['Window Size'], w['Hit Count']) for w in matched_windows}
        size_hit_str = ", ".join(f"(W{w},H{h})" for w, h in sorted(configs))

        # Find the window with max hit/window ratio
        best_ratio = 0.0
        best_window_seq = ""
        for w_info in matched_windows:
            ratio = w_info['Hit Count'] / w_info['Window Size']
            if ratio > best_ratio:
                best_ratio = ratio
                best_window_seq = w_info['Window Sequence']

        max_rate_window_value = f"{best_ratio:.4f}: {best_window_seq}"

        row = [
            gene_data['Gene Name'],
            gene_data['Sequence'],
            window_str,
            size_hit_str,
            max_rate_window_value
        ]
        ws_unique.append(row)

    wb.save(excel_file)
    print(f"Task '{task_words}' results saved to: {excel_file}")
    print(f"   - Number of unique sequences: {len(unique_sequences_entries)}")


def main():
    parser = argparse.ArgumentParser(description='Pattern matching analysis of gene sequences (support wildcard .)')
    parser.add_argument('tasks', type=str, help='Task list, e.g.: ATCG,GCTA or D-E,F-G')
    parser.add_argument('--input', type=str, default='./Database/normal.fasta', help='Input FASTA file path (support .gz)')
    parser.add_argument('--start', type=int, default=20,help='Start window size (optional)')
    parser.add_argument('--end', type=int, default=100, help='End window size (default 100)')
    parser.add_argument('--output', type=str, default='./Result', help='Output directory (default ./Result)')
    parser.add_argument('--point', nargs='+', type=int, metavar=('HIT', 'WINDOW'), help='Specify points, e.g. 2 20 or 2 20 4 80')
    parser.add_argument('--mode', choices=['formula', 'fix'], default='formula', help='Mode: formula (default) or fix')
    parser.add_argument('--outputtype', choices=['all', 'simple'], default='simple',
                        help='Output type: all (include all hits) or simple (only unique sequences)')

    args = parser.parse_args()

    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")

    sequences = load_sequences(args.input)
    print(f"Loaded {len(sequences)} sequences from {args.input}.")

    os.makedirs(args.output, exist_ok=True)

    if args.mode == 'formula':
        if args.point is None:
            hit1, window1 = 18, 20
            hit2, window2 = 30, 50
            points = (hit1, window1, hit2, window2)
            k, b = calculate_k_b(hit1, window1, hit2, window2)
            print(f"Warning: --point not provided, using default: k={k:.6f}, b={b:.6f} "
                  f"(based on points (hit={hit1}->window={window1}) and (hit={hit2}->window={window2}))")
        elif len(args.point) == 2:
            h1, w1 = args.point
            points = (h1, w1)
            k, b = calculate_k_b(h1, w1)
            print(f"Single point mode: hit={h1} -> window={w1}, k={k:.6f}, b={b:.6f}")
        elif len(args.point) == 4:
            h1, w1, h2, w2 = args.point
            points = (h1, w1, h2, w2)
            k, b = calculate_k_b(h1, w1, h2, w2)
            print(f"Two point mode: ({h1},{w1}), ({h2},{w2}) -> k={k:.6f}, b={b:.6f}")
        else:
            raise ValueError("--point argument count error, should be 2 or 4 integers")

        start_window = args.start if args.start is not None else 10

    else:
        if args.point is None:
            raise ValueError("In --mode fix, --point parameter must be provided")
        points = tuple(args.point)
        k, b = None, None
        start_window = args.start if args.start is not None else 1

    end_window = args.end
    if start_window > end_window:
        raise ValueError(f"Start window ({start_window}) cannot be greater than end window ({end_window})")

    print(f"Window range: {start_window} ~ {end_window}")

    for task in args.tasks.split(','):
        task = task.strip()
        if not task:
            continue
        handle_task(
            sequences=sequences,
            task_words=task,
            start_window=start_window,
            end_window=end_window,
            mode=args.mode,
            points=points,
            input_filename=args.input,
            output_dir=args.output,
            outputtype=args.outputtype,
            k_val=k,
            b_val=b
        )


if __name__ == "__main__":
    main()