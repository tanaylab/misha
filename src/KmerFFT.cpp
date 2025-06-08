#include "KmerFFT.h"
#include "GenomeSeqFetch.h"
#include <algorithm>
#include <numeric>
#include <limits>
#include <cctype>

// KISS FFT Library
// Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
// Licensed under BSD-3-Clause - https://github.com/mborgerding/kissfft

KmerFFT::KmerFFT(const std::string& kmer, const std::string& groot, Mode mode, 
                 bool extend, double freq, WindowType window)
    : GenomeSeqScorer(groot, extend, 0), m_kmer(kmer), m_mode(mode), 
      m_freq(freq), m_window(window) {
    m_kmer_len = kmer.length();
    // Convert kmer to uppercase
    std::transform(m_kmer.begin(), m_kmer.end(), m_kmer.begin(), ::toupper);
}

float KmerFFT::score_interval(const GInterval& interval, const GenomeChromKey& chromkey) {
    // Calculate the appropriate interval based on extension setting
    GInterval expanded_interval = calculate_expanded_interval(interval, chromkey, m_kmer_len);
    
    if (expanded_interval.start >= expanded_interval.end && (interval.end - interval.start) > 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }
    
    std::vector<char> seq_vec;
    if (expanded_interval.end > expanded_interval.start) {
        m_seqfetch.read_interval(expanded_interval, chromkey, seq_vec);
    }
    
    if (seq_vec.empty() && (interval.end - interval.start) > 0) {
        return std::numeric_limits<float>::quiet_NaN();
    }
    
    std::string seq(seq_vec.begin(), seq_vec.end());
    
    size_t signal_len = interval.end - interval.start;
    size_t seq_offset = m_extend ? interval.start - expanded_interval.start : 0;

    // Get the kmer occurrence signal
    std::vector<double> signal = get_kmer_signal(seq, signal_len, seq_offset);
    
    if (signal.size() < 2) {
        if (m_mode == POWER_AT_FREQ) {
            return 0;
        }
        return std::numeric_limits<float>::quiet_NaN();
    }
    
    // Apply window function
    apply_window(signal);
    
    // Compute FFT
    std::vector<std::complex<double>> fft_result;
    compute_fft(signal, fft_result);
    
    // Compute power spectrum
    std::vector<double> power_spectrum(fft_result.size() / 2 + 1);
    for (size_t i = 0; i < power_spectrum.size(); ++i) {
        power_spectrum[i] = std::norm(fft_result[i]) / signal.size();
    }
    
    switch (m_mode) {
        case POWER_AT_FREQ: {
            // Convert frequency to bin index
            double bin_index = m_freq * signal.size();
            int bin = static_cast<int>(std::round(bin_index));
            
            if (bin < 0 || bin >= (int)power_spectrum.size()) {
                return 0.0;
            }
            
            return static_cast<float>(power_spectrum[bin]);
        }
        
        case PEAK_FREQ: {
            // Find the peak frequency (skip DC component)
            if (power_spectrum.size() <= 1) {
                return 0.0;
            }
            auto max_it = std::max_element(power_spectrum.begin() + 1, power_spectrum.end());
            size_t max_idx = std::distance(power_spectrum.begin(), max_it);
            
            // Convert bin index to frequency (cycles per base)
            return static_cast<float>(static_cast<double>(max_idx) / signal.size());
        }
        
        case PEAK_POWER: {
            // Find the peak power (skip DC component)
            if (power_spectrum.size() <= 1) {
                return 0.0;
            }
            auto max_it = std::max_element(power_spectrum.begin() + 1, power_spectrum.end());
            return static_cast<float>(*max_it);
        }
    }
    
    return std::numeric_limits<float>::quiet_NaN();
}

std::vector<double> KmerFFT::get_kmer_signal(const std::string& sequence, size_t signal_len, size_t seq_offset) {
    std::vector<double> signal(signal_len, 0.0);

    if (m_kmer_len == 0 || signal_len == 0) {
        return signal;
    }

    for (size_t i = 0; i < signal_len; ++i) {
        size_t seq_idx = i + seq_offset;
        if (seq_idx + m_kmer_len > sequence.length()) {
            break;
        }
        bool match = true;
        for (size_t j = 0; j < m_kmer_len; ++j) {
            char seq_char = std::toupper(sequence[seq_idx + j]);
            if (seq_char != m_kmer[j]) {
                match = false;
                break;
            }
        }
        if (match) {
            signal[i] = 1.0;
        }
    }
    
    return signal;
}

void KmerFFT::apply_window(std::vector<double>& signal) {
    size_t n = signal.size();
    
    switch (m_window) {
        case WINDOW_NONE:
            break;
            
        case WINDOW_HANN:
            for (size_t i = 0; i < n; ++i) {
                double w = 0.5 * (1 - std::cos(2 * M_PI * i / (n - 1)));
                signal[i] *= w;
            }
            break;
            
        case WINDOW_HAMMING:
            for (size_t i = 0; i < n; ++i) {
                double w = 0.54 - 0.46 * std::cos(2 * M_PI * i / (n - 1));
                signal[i] *= w;
            }
            break;
            
        case WINDOW_BLACKMAN:
            for (size_t i = 0; i < n; ++i) {
                double w = 0.42 - 0.5 * std::cos(2 * M_PI * i / (n - 1)) + 
                          0.08 * std::cos(4 * M_PI * i / (n - 1));
                signal[i] *= w;
            }
            break;
    }
}

// KISS FFT implementation
void KmerFFT::compute_fft(const std::vector<double>& signal, 
                          std::vector<std::complex<double>>& fft_result) {
    int n = signal.size();
    fft_result.resize(n);
    
    if (n <= 1) {
        if (n == 1) {
            fft_result[0] = std::complex<double>(signal[0], 0.0);
        }
        return;
    }
    
    // Allocate KISS FFT configuration
    kiss_fft_cfg cfg = kiss_fft_alloc(n, 0, NULL, NULL); // 0 = forward FFT
    if (!cfg) {
        // Fallback to simple DFT if allocation fails
        for (int k = 0; k < n; ++k) {
            std::complex<double> sum(0, 0);
            for (int t = 0; t < n; ++t) {
                double angle = -2 * M_PI * k * t / n;
                sum += signal[t] * std::complex<double>(std::cos(angle), std::sin(angle));
            }
            fft_result[k] = sum;
        }
        return;
    }
    
    // Prepare input data for KISS FFT (convert to kiss_fft_cpx format)
    std::vector<kiss_fft_cpx> cx_in(n), cx_out(n);
    for (int i = 0; i < n; ++i) {
        cx_in[i].r = (kiss_fft_scalar)signal[i];
        cx_in[i].i = 0.0;
    }
    
    // Perform FFT
    kiss_fft(cfg, cx_in.data(), cx_out.data());
    
    // Convert results back to std::complex<double>
    for (int i = 0; i < n; ++i) {
        fft_result[i] = std::complex<double>(cx_out[i].r, cx_out[i].i);
    }
    
    // Clean up
    kiss_fft_free(cfg);
} 