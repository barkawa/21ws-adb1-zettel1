import numpy as np
import matplotlib.pyplot as plt
import re

L = 30
START_CODON_IDX = 100
PSEUDOCOUNT = 1
BACKGROUND = 0.25
FILEPATH = '../data/TIS-Ecoli.txt'



def base_to_pwm_row(base):
    ''' Maps each base to a numerical index, representing a matrix row '''
    return {'A': 0, 'C': 1, 'T': 2, 'G': 3}[base]


def count_start_codons(samples):
    n = 0
    for sample in samples:
        n += sum(1 for _ in re.finditer(r'(?=(ATG|GTG|TTG))', sample))
    return n


def calc_pfm(samples):
    ''' Calculate a Position Frequency Matrix (PFM)
    for a region of L bases before the start codon,
    using a set of samples '''

    pfm = np.zeros((4, L))

    for sample in samples:
        for i in range(0, L):
            base = sample[START_CODON_IDX - L + i]
            pfm[base_to_pwm_row(base), i] += 1

    return pfm


class GeneStartCandidate:
    ''' A possible gene start, identified by looking for ATG/GTG/TTG in a sequence
    and scored with a PWM '''
    def __init__(self, sample_index, start_index, score):
        self.sample_index = sample_index
        self.start_index = start_index
        self.score = score


def find_and_score_gene_start_candidates(samples, pwm):
    ''' Look for possible gene starts, and score them using a PWM.
    only using candidates, for which all positions of the PWM are observable '''

    gene_start_candidates = []

    for sample_index, sample in enumerate(samples):
        for candidate in re.finditer(r'(?=(ATG|GTG|TTG))', sample):
            if candidate.start() - L >= 0:
                score = 0
                for i in range(0, L):
                        base = sample[candidate.start() - L + i]
                        score += pwm[base_to_pwm_row(base), i]
                gene_start_candidates.append(GeneStartCandidate(sample_index, candidate.start(), score))
    
    return gene_start_candidates


def estimate_score_threshold(candidates, sensitivity_percent):
    ''' Estimate a score threshold for a desired sensitivity,
    by looking at the distribution of scores of real positive cases '''

    real_positives = [c.score for c in candidates if c.start_index == 100]
    return np.percentile(real_positives, sensitivity_percent)


def calc_pwm(samples):
    # Calculate the Position Frequency Matrix (PFM) for the given sample set
    # including the pseudocount
    pfm = calc_pfm(samples)
    pfm += PSEUDOCOUNT

    # Estimate the Position Probability Matrix (PPM)
    # by dividing the elements of the PFM by the sample count
    # (also taking the pseudocount into consideration)
    ppm = pfm / (len(samples) + PSEUDOCOUNT * 4)

    # Estimate the Position Weight Matrix (PWM) by
    # calculating the difference of log-likelihoods
    # between the PPM and a background model (0.25, ..., 0.25)
    pwm = np.log2(ppm) - np.log2(BACKGROUND)

    return pwm


if __name__ == "__main__":

    # Open the file and create a list of DNA samples from it
    with open(FILEPATH, "r") as file:
        samples = [line.rstrip() for line in file.readlines()]
    
    # Count all start codons in the dataset
    print("Anzahl der Startcodons: ", count_start_codons(samples))

    print(" --- Trainingsdaten = Validierungsdaten --- ")

    pwm = calc_pwm(samples)

    # Plot the PWM
    plt.figure(1, figsize=(8,3))
    plt.imshow(pwm, cmap='viridis', interpolation='nearest')
    plt.colorbar()
    plt.yticks(np.arange(0,4), ['A','C','T','G'])
    plt.savefig('../out/plot1.pdf')

    # Look for start codons (= gene start candidates) and score them
    candidates = find_and_score_gene_start_candidates(samples, pwm)

    # Estimate threshold for a 50% target sensitivity
    threshold = estimate_score_threshold(candidates, 50)
    print("Score-Threshold für eine Sensitivität von 50%: ", threshold)

    # Get all matches above the threshold
    positives = [c for c in candidates if c.score > threshold]

    # Discriminate true and false positives
    true_positives = [c for c in positives if c.start_index == 100]
    false_positives = [c for c in positives if c.start_index != 100]
    print("Echt/falsch positive Kandidaten:", len(true_positives), "/", len(false_positives))


    #---------------------------------------------------#
    print(" --- Trainingsd. von Validierungsd. getrennt--- ")
    # Split samples into training and validation sets
    training = samples[:400]
    validation = samples[400:]

    pwm = calc_pwm(training)
    training_candidates = find_and_score_gene_start_candidates(training, pwm)
    threshold = estimate_score_threshold(training_candidates, 50)
    print("Threshold: ", threshold)

    validation_candidates = find_and_score_gene_start_candidates(validation, pwm)

    # Get all matches above the threshold
    positives = [c for c in validation_candidates if c.score > threshold]

    # Discriminate true and false positives
    true_positives = [c for c in positives if c.start_index == 100]
    false_positives = [c for c in positives if c.start_index != 100]
    print("Echt/falsch positive Kandidaten:", len(true_positives), "/", len(false_positives))

    # ROC Plot
    p = len(validation)
    n = count_start_codons(validation) - p
    threshold_range = np.flip(np.arange(-100., 100., 0.1))
    tprs = []
    fprs = []
    for t in threshold_range:
        positives = [c for c in validation_candidates if c.score > t]
        tp = sum(1 for _ in [c for c in positives if c.start_index == 100])
        fp = sum(1 for _ in [c for c in positives if c.start_index != 100])
        tprs.append(tp / p)
        fprs.append(fp / n)

    plt.figure(2)
    plt.plot(fprs, tprs)
    plt.plot([0, 1], 'k:')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.gca().set_aspect('equal')
    plt.savefig('../out/roc.pdf')

    print("Area under Curve: ", np.trapz(tprs, fprs))
