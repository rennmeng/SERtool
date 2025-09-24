use pyo3::prelude::*;
use pyo3::types::PyDict;
use rayon::prelude::*;
use std::collections::HashMap;

// 引入缺失的异常类型
use pyo3::exceptions::PyValueError;

#[pyclass]
#[derive(Clone, Debug)]
pub struct HitInfo {
    #[pyo3(get, set)]
    pub gene_name: String,
    #[pyo3(get, set)]
    pub window_start: usize,
    #[pyo3(get, set)]
    pub window_end: usize,
    #[pyo3(get, set)]
    pub window_sequence: String,
    #[pyo3(get, set)]
    pub window_size: usize,
    #[pyo3(get, set)]
    pub hit_count: u32,
}

fn matches_pattern(pattern: &str, text: &str) -> bool {
    if pattern.len() != text.len() {
        return false;
    }
    pattern
        .chars()
        .zip(text.chars())
        .all(|(p, t)| p == '.' || p == t)
}

fn count_hits_in_window(window: &str, words: &[&str]) -> u32 {
    let mut total = 0;
    for &word in words {
        if word.is_empty() {
            continue;
        }
        if word.contains('.') {
            let w_len = word.len();
            if w_len <= window.len() {
                for i in 0..=(window.len() - w_len) {
                    if matches_pattern(word, &window[i..i + w_len]) {
                        total += 1;
                    }
                }
            }
        } else {
            let mut start = 0;
            while let Some(pos) = window[start..].find(word) {
                total += 1;
                start += pos + 1;
            }
        }
    }
    total
}

#[pyfunction]
fn scan_sequences_batch(
    _py: Python,  // 显式传入但未使用，加 _ 前缀避免警告
    sequences: &Bound<'_, PyDict>,
    words: Vec<String>,
    window_size: usize,
    min_hit: u32,
) -> PyResult<Vec<HitInfo>> {
    // Step 1: 将 PyDict 复制到 Rust HashMap，脱离 GIL
    let seq_map: HashMap<String, String> = sequences
        .iter()
        .map(|(k, v)| {
            let key: String = k.extract().map_err(|_| {
                PyErr::new::<PyValueError, _>("Dictionary key must be string")
            })?;
            let value: String = v.extract().map_err(|_| {
                PyErr::new::<PyValueError, _>("Sequence value must be string")
            })?;
            Ok((key, value))
        })
        .collect::<PyResult<_>>()?;

    // Step 2: 转为 Vec 以便并行处理
    let seq_vec: Vec<(String, String)> = seq_map.into_iter().collect();
    let word_slices: Vec<&str> = words.iter().map(|s| s.as_str()).collect();

    // Step 3: 并行扫描（完全在 Rust 中，无 GIL）
    let results: Vec<HitInfo> = seq_vec
        .into_par_iter()
        .flat_map(|(gene_name, sequence)| {
            let mut gene_hits = Vec::new();

            if sequence.len() < window_size {
                return gene_hits;
            }

            for i in 0..=(sequence.len() - window_size) {
                let window = &sequence[i..i + window_size];
                let hit_count = count_hits_in_window(window, &word_slices);
                if hit_count >= min_hit {
                    gene_hits.push(HitInfo {
                        gene_name: gene_name.clone(),
                        window_start: i,
                        window_end: i + window_size,
                        window_sequence: window.to_string(),
                        window_size,
                        hit_count,
                    });
                }
            }

            gene_hits
        })
        .collect();

    Ok(results)
}

#[pyfunction]
fn calculate_k_b(
    hit1: u32,
    window1: usize,
    hit2: Option<u32>,
    window2: Option<usize>,
) -> PyResult<(f64, f64)> {
    if hit1 == 0 {
        return Err(PyErr::new::<PyValueError, _>("Hit1 must be > 0"));
    }

    let (k, b) = match (hit2, window2) {
        (Some(h2), Some(w2)) => {
            if h2 == 0 || h2 == hit1 {
                return Err(PyErr::new::<PyValueError, _>("Invalid second point"));
            }
            let k = (w2 as f64 - window1 as f64) / (h2 as f64 - hit1 as f64);
            let b = window1 as f64 - k * hit1 as f64;
            (k, b)
        }
        (None, None) => (window1 as f64 / hit1 as f64, 0.0),
        _ => {
            return Err(PyErr::new::<PyValueError, _>(
                "Both hit2 and window2 must be provided",
            ));
        }
    };

    Ok((k, b))
}

// ✅ 必须添加：模块入口
#[pymodule]
fn ser_tool(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction_bound!(scan_sequences_batch, m)?)?;
    m.add_function(wrap_pyfunction_bound!(calculate_k_b, m)?)?;
    m.add_class::<HitInfo>()?;
    Ok(())
}