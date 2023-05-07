//! Rohde & Schwarz Engineering Competition 2023
//!
//! This is the code to speed up. Enjoy!

#pragma once

#include "ec2023/ec2023.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

static constexpr float OVERLAP_RATIO = 0.75;
static constexpr size_t WINDOW_SIZE = 1024;
static constexpr size_t SIZE_SPECTRUM = (WINDOW_SIZE / 2) + 1;
static const ec::Float PI = 3.14159265358979323846f;

/*
void compute_fourier_transform(const std::vector<ec::Float> &input, std::vector<ec::Float> &outputReal,
                               std::vector<ec::Float> &outputImag);// */

void compute_fourier_transform(const std::vector<ec::Float> &input, std::vector<ec::Float> &outputReal,
                               std::vector<ec::Float> &outputImag) {

    //const ec::Float PI = 3.14159265358979323846f;

    size_t inputSize = input.size();

    outputReal.clear();
    outputReal.resize(inputSize, 0.0f);
    outputImag.clear();
    outputImag.resize(inputSize, 0.0f);

    //TODO: Wieder Nestedloops
    for (size_t I = 0; I < inputSize; ++I) {
        const ec::Float angleTerm = (-2.0f * PI) * ec::Float(I) * (1.0f / ec::Float(inputSize));

        for (size_t J = 0; J < inputSize; ++J) {
            //const ec::Float angleTerm = (-2.0f * PI) * ec::Float(I) * (1.0f / ec::Float(inputSize));

            outputReal[I] += input[J] * ec_cos(angleTerm * J);
            outputImag[I] += input[J] * ec_sin(angleTerm * J);
        }
    }
}

std::vector<ec::Float> process_signal(const std::vector<ec::Float> &inputSignal) {
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

    //TODO: Optimierungspotential -NestedLoops-
    for (size_t J = 0; J < numWins; J++) {
        for (size_t I = 0; I < WINDOW_SIZE; I++) {

            //TODO: Find a better way to declare this variable
            ec::Float blackmanWinCoef = 0.42f - 0.5f * ec_cos(ec::Float(I) * 2.0f * PI / (WINDOW_SIZE - 1));
            blackmanWinCoef += 0.08f * ec_cos(ec::Float(I) * 4.0f * PI / (WINDOW_SIZE - 1));

            signalWindow[I] = inputSignal[I + idxStartWin] * blackmanWinCoef;
        }

        compute_fourier_transform(signalWindow, signalFreqReal, signalFreqImag);

        for (size_t I = 0; I < SIZE_SPECTRUM; I++) {

            ec::Float freqVal = signalFreqReal[I] * signalFreqReal[I] + signalFreqImag[I] * signalFreqImag[I];
            freqVal = ec_sqrt(freqVal);
            freqVal = freqVal / ec::Float(WINDOW_SIZE);

            if (I > 0 && I < SIZE_SPECTRUM - 1)
                freqVal = freqVal * 2.0f;

            freqVal = freqVal * freqVal;

            freqVal = 10.0f * ec_log10(1000.0f * freqVal);

            outputSpectrum[I] = ec_max(outputSpectrum[I], freqVal);
        }

        idxStartWin += stepBetweenWins;

    } //TODO: Endet hier....

    return outputSpectrum;
}// */

