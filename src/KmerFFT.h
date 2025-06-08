#ifndef KMERFFT_H_
#define KMERFFT_H_

#include <string>
#include <vector>
#include <complex>
#include <cmath>
#include <memory>
#include "GInterval.h"
#include "GenomeChromKey.h"
#include "GenomeSeqScorer.h"
#include "kiss_fft.h"



class KmerFFT : public GenomeSeqScorer {
public:
    enum Mode {
        POWER_AT_FREQ,
        PEAK_FREQ,
        PEAK_POWER
    };

    enum WindowType {
        WINDOW_NONE,
        WINDOW_HANN,
        WINDOW_HAMMING,
        WINDOW_BLACKMAN
    };

    KmerFFT(const std::string& kmer, const std::string& groot, Mode mode, 
            bool extend = true, double freq = 0.0, WindowType window = WINDOW_HANN);

    float score_interval(const GInterval& interval, const GenomeChromKey& chromkey) override;

private:
    std::string m_kmer;
    Mode m_mode;
    double m_freq;
    WindowType m_window;
    size_t m_kmer_len;

    void compute_fft(const std::vector<double>& signal, std::vector<std::complex<double>>& fft_result);
    void apply_window(std::vector<double>& signal);
    std::vector<double> get_kmer_signal(const std::string& sequence, size_t signal_len, size_t seq_offset);
};

#endif /* KMERFFT_H_ */ 