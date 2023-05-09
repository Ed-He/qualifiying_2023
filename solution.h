//! Rohde & Schwarz Engineering Competition 2023
//!
//! This is the code to speed up. Enjoy!

#pragma once

#include "ec2023/ec2023.h"
#include <iostream>      
#include <iomanip>
#include <vector>

static constexpr float OVERLAP_RATIO = 0.75;
static constexpr size_t WINDOW_SIZE = 1024;
static constexpr size_t SIZE_SPECTRUM = (WINDOW_SIZE / 2) + 1;
static const ec::Float PI = 3.14159265358979323846f;
static const ec::Float WINDOW1 = PI / (WINDOW_SIZE - 1) * 2.0f;
static const ec::Float WINDOW2 = PI / (WINDOW_SIZE - 1) * 4.0f;

void compute_fourier_transform(const std::vector<ec::Float>& input, std::vector<ec::Float>& outputReal, std::vector<ec::Float>& outputImag);
size_t bit_reverse(size_t n, size_t bits);


void compute_fourier_transform(const std::vector<ec::Float>& input, std::vector<ec::Float>& outputReal, std::vector<ec::Float>& outputImag)
{
    const ec::Float PI = 3.14159265358979323846f;

    size_t inputSize = input.size();

    outputReal.clear();
    outputReal.resize(inputSize, 0.0f);
    outputImag.clear();
    outputImag.resize(inputSize, 0.0f);

    // �berpr�fe, ob die Eingabegr��e eine Potenz von 2 ist
    if ((inputSize & (inputSize - 1)) != 0) {
        std::cerr << "Error: Input size must be a power of 2" << std::endl;
        return;
    }

    // Berechne die FFT mit der Cooley-Tukey-Methode
    for (size_t i = 0; i < inputSize; ++i) {
        size_t j = bit_reverse(i, log2(inputSize)); // Bit-Umkehrung

        outputReal[j] = input[i]; // Initialisiere mit der Eingabe
        outputImag[j] = 0.0;
    }

    for (size_t size = 2; size <= inputSize; size *= 2) { // Iteriere durch die verschiedenen Stufen der FFT
        size_t half_size = size / 2;
        ec::Float angle = 2 * PI / size;
        ec::Float w_real = ec::ec_cos(angle);
        ec::Float w_imag = ec::ec_sin(angle);

        for (size_t i = 0; i < inputSize; i += size) { // Iteriere durch die Bl�cke der Gr��e "size"
            ec::Float wr = 1.0f;
            ec::Float wi = 0.0f;

            for (size_t j = 0; j < half_size; ++j) { // Iteriere durch die Elemente innerhalb jedes Blocks
                ec::Float re = outputReal[i + j + half_size] * wr - outputImag[i + j + half_size] * wi;
                ec::Float im = outputReal[i + j + half_size] * wi + outputImag[i + j + half_size] * wr;

                outputReal[i + j + half_size] = outputReal[i + j] - re;
                outputImag[i + j + half_size] = outputImag[i + j] - im;
                outputReal[i + j] += re;
                outputImag[i + j] += im;

                ec::Float t = wr;
                wr = wr * w_real - wi * w_imag;
                wi = t * w_imag + wi * w_real;
            }
        }
    }
}

// Funktion zur Berechnung der Bit-Umkehrung
size_t bit_reverse(size_t n, size_t bits) {
    size_t reversed_n = n;
    size_t count = bits - 1;

    n >>= 1;
    while (n) {
        reversed_n <<= 1;
        reversed_n |= n & 1;
        n >>= 1;
        --count;
    }

    return (reversed_n << count) & ((1 << bits) - 1);
}

/*std::vector<ec::Float> process_signal(const std::vector<ec::Float>& inputSignal)
{  //kfadhsjfnab
  const size_t numSamples = inputSignal.size();
  const size_t sizeSpectrum = (WINDOW_SIZE / 2) + 1;
  const size_t stepBetweenWins = static_cast<size_t>(ceil(WINDOW_SIZE * (1 - OVERLAP_RATIO)));
  const size_t numWins = (numSamples - WINDOW_SIZE) / stepBetweenWins + 1;
  const ec::Float PI = 3.14159265358979323846f;

  std::vector<ec::Float> signalWindow(WINDOW_SIZE);
  std::vector<ec::Float> signalFreqReal(WINDOW_SIZE);
  std::vector<ec::Float> signalFreqImag(WINDOW_SIZE);
  std::vector<ec::Float> spectrumWindow(sizeSpectrum);
  std::vector<ec::Float> outputSpectrum(sizeSpectrum, std::numeric_limits<float>::lowest());

  size_t idxStartWin = 0;

  for (size_t J = 0; J < numWins; J++)
  {
    for (size_t I = 0; I < WINDOW_SIZE; I++)
    {
      ec::Float blackmanWinCoef = 0.42f - 0.5f * ec_cos(ec::Float(I) * 2.0f * PI / (WINDOW_SIZE - 1));                                         
      blackmanWinCoef = blackmanWinCoef + 0.08f * ec_cos(ec::Float(I) * 4.0f * PI / (WINDOW_SIZE - 1));

      signalWindow[I] = inputSignal[I + idxStartWin] * blackmanWinCoef;
    }

    compute_fourier_transform(signalWindow, signalFreqReal, signalFreqImag);

    for (size_t I = 0; I < sizeSpectrum; I++)
    {
      ec::Float freqVal = signalFreqReal[I] * signalFreqReal[I] + signalFreqImag[I] * signalFreqImag[I];
      freqVal = ec_sqrt(freqVal);
      freqVal = freqVal / ec::Float(WINDOW_SIZE);

      if (I > 0 && I < sizeSpectrum - 1) freqVal = freqVal * 2.0f;

      freqVal = freqVal * freqVal;

      freqVal = 10.0f * ec_log10(1000.0f * freqVal);

      outputSpectrum[I] = ec_max(outputSpectrum[I], freqVal);
    }

    idxStartWin += stepBetweenWins;

  }

  return outputSpectrum;
}*/
std::vector<ec::Float> process_signal(const std::vector<ec::Float>& inputSignal) {
    const size_t numSamples = inputSignal.size();
    //const size_t sizeSpectrum = (WINDOW_SIZE / 2) + 1;
    const size_t stepBetweenWins = static_cast<size_t>(std::ceil(WINDOW_SIZE * (1 - OVERLAP_RATIO)));
    const size_t numWins = (numSamples - WINDOW_SIZE) / stepBetweenWins + 1;
    //const ec::Float PI = 3.14159265358979323846f;

    std::vector<ec::Float> signalWindow(WINDOW_SIZE);
    std::vector<ec::Float> signalFreqReal(WINDOW_SIZE);
    std::vector<ec::Float> signalFreqImag(WINDOW_SIZE);
    //std::vector<ec::Float> spectrumWindow(SIZE_SPECTRUM);
    std::vector<ec::Float> outputSpectrum(SIZE_SPECTRUM, std::numeric_limits<float>::lowest());

    size_t idxStartWin = 0;

    ec::Float blackmanWinCoef;
    ec::Float freqVal;
    ec::Float freqWindow = ec::Float(WINDOW_SIZE);

    //TODO: Optimierungspotential -NestedLoops-
    for (size_t J = 0; J < numWins; J++) {
        for (size_t I = 0; I < WINDOW_SIZE; I++) {

            //TODO: Find a better way to declare this variable
            blackmanWinCoef = 0.42f - (0.5f * ec_cos(ec::Float(I) * WINDOW1)) + (0.08f * ec_cos(ec::Float(I) * WINDOW2));
            //blackmanWinCoef = 0.42f - 0.5f * ec_cos(ec::Float(I) * 2.0f * PI / (WINDOW_SIZE - 1));
            //blackmanWinCoef += 0.08f * ec_cos(ec::Float(I) * 4.0f * PI / (WINDOW_SIZE - 1));

            signalWindow[I] = inputSignal[I + idxStartWin] * blackmanWinCoef;
        }

        idxStartWin += stepBetweenWins;

        compute_fourier_transform(signalWindow, signalFreqReal, signalFreqImag);

        freqVal = ec_sqrt(signalFreqReal[0] * signalFreqReal[0] + signalFreqImag[0] * signalFreqImag[0]) / freqWindow;
        freqVal = 30.0f + 10.0f * ec_log10(freqVal * freqVal);
        outputSpectrum[0] = ec_max(outputSpectrum[0], freqVal);

        for (size_t I = 1; I < SIZE_SPECTRUM - 1; I++) {

            freqVal = ec_sqrt(signalFreqReal[I] * signalFreqReal[I] + signalFreqImag[I] * signalFreqImag[I]) / freqWindow;
            freqVal = 30.0f + 10.0f * ec_log10(freqVal * freqVal * 4.0f);
            outputSpectrum[I] = ec_max(outputSpectrum[I], freqVal);
        }

        size_t L = SIZE_SPECTRUM - 1;
        freqVal = ec_sqrt(signalFreqReal[L] * signalFreqReal[L] + signalFreqImag[L] * signalFreqImag[L]) / freqWindow;
        freqVal = 30.0f + 10.0f * ec_log10(freqVal * freqVal);
        outputSpectrum[L] = ec_max(outputSpectrum[L], freqVal);

        //idxStartWin += stepBetweenWins;

    }//TODO: Endet hier....

    return outputSpectrum;
}// */





