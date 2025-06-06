#include "KmerFFT.h"
#include "GenomeSeqFetch.h"
#include <algorithm>
#include <numeric>
#include <limits>
#include <cctype>

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
    
    if (expanded_interval.start >= expanded_interval.end) {
        return std::numeric_limits<float>::quiet_NaN();
    }
    
    std::vector<char> seq_vec;
    m_seqfetch.read_interval(expanded_interval, chromkey, seq_vec);
    
    if (seq_vec.empty()) {
        return std::numeric_limits<float>::quiet_NaN();
    }
    
    std::string seq(seq_vec.begin(), seq_vec.end());
    
    // Get the kmer occurrence signal
    std::vector<double> signal = get_kmer_signal(seq);
    
    if (signal.size() < 2) {
        return std::numeric_limits<float>::quiet_NaN();
    }
    
    // Apply window function
    apply_window(signal);
    
    // Compute FFT
    std::vector<std::complex<double>> fft_result;
    compute_fft(signal, fft_result);
    
    // Compute power spectrum
    std::vector<double> power_spectrum(fft_result.size() / 2);
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
            auto max_it = std::max_element(power_spectrum.begin() + 1, power_spectrum.end());
            size_t max_idx = std::distance(power_spectrum.begin(), max_it);
            
            // Convert bin index to frequency (cycles per base)
            return static_cast<float>(static_cast<double>(max_idx) / signal.size());
        }
        
        case PEAK_POWER: {
            // Find the peak power (skip DC component)
            auto max_it = std::max_element(power_spectrum.begin() + 1, power_spectrum.end());
            return static_cast<float>(*max_it);
        }
    }
    
    return std::numeric_limits<float>::quiet_NaN();
}

std::vector<double> KmerFFT::get_kmer_signal(const std::string& sequence) {
    std::vector<double> signal(sequence.length(), 0.0);
    
    for (size_t i = 0; i <= sequence.length() - m_kmer_len; ++i) {
        bool match = true;
        for (size_t j = 0; j < m_kmer_len; ++j) {
            char seq_char = std::toupper(sequence[i + j]);
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

// Simple DFT implementation
void KmerFFT::compute_fft(const std::vector<double>& signal, 
                          std::vector<std::complex<double>>& fft_result) {
    size_t n = signal.size();
    fft_result.resize(n);
    
    for (size_t k = 0; k < n; ++k) {
        std::complex<double> sum(0, 0);
        for (size_t t = 0; t < n; ++t) {
            double angle = -2 * M_PI * k * t / n;
            sum += signal[t] * std::complex<double>(std::cos(angle), std::sin(angle));
        }
        fft_result[k] = sum;
    }
} 