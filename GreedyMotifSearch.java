import java.util.ArrayList;
import java.util.List;

import edu.princeton.cs.algs4.In;

/**
 * 
 * A solution to a programming assignment for the Bioinformatics Algorithms
 * (Part 1) on Coursera. The associated textbook is Bioinformatics Algorithms:
 * An Active-Learning Approach by Phillip Compeau & Pavel Pevzner. The course is
 * run on Coursera and the assignments and textbook are hosted on Stepic
 * 
 * Problem Title: Greedy Motif Search 
 * Assignment #: 03 
 * Problem ID: D 
 * URL:https://beta.stepic.org/Bioinformatics-Algorithms-2/Greedy-Motif-Search-159/#step-5
 * 
 */
public class GreedyMotifSearch {

    private List<String> bestMotifs;

    private double bestScore;

    private double[][] currentProfile;

    private boolean applyLaplacian;

    public GreedyMotifSearch(List<String> dNaList, int k, int t, boolean applyLaplacian) {
        this.applyLaplacian = applyLaplacian;
        bestScore = Double.MAX_VALUE;
        String firstDna = dNaList.get(0);
        for (int i = 0; i <= firstDna.length() - k; i++) {
            List<String> currentMotifs = new ArrayList<String>();
            currentMotifs.add(firstDna.substring(i, i + k));

            currentProfile = profile(currentMotifs);
            for (int j = 1; j < t; j++) {
                currentMotifs.add(Bioinformatics.profileMostProbable(dNaList.get(j), k, currentProfile));
                currentProfile = profile(currentMotifs);
            }

            double currentScore = score(currentMotifs);

            if (currentScore < bestScore) {
                bestScore = currentScore;
                bestMotifs = new ArrayList<String>(currentMotifs);
            }
        }
    }

    private List<String> getBestMotifs() {
        return bestMotifs;
    }

    private int score(List<String> motifs) {
        double[][] profileMatrix = profile(motifs);
        int score = 0;

        for (int i = 0; i < motifs.get(0).length(); i++) {
            char maxChar = getMaxOccurinChar(profileMatrix, i);
            for (int j = 0; j < 4; j++) {
                if (j != getNum(maxChar)) {
                    score += profileMatrix[j][i];
                }
            }
        }

        return score;
    }

    private char getMaxOccurinChar(double[][] profileMatrix, int column) {
        double maxNum = Double.MIN_VALUE;
        int maxCharId = -1;

        for (int i = 0; i < 4; i++) {
            if (profileMatrix[i][column] > maxNum) {
                maxNum = profileMatrix[i][column];
                maxCharId = i;
            }
        }

        return getProtien(maxCharId);
    }

    private double[][] profile(List<String> motifs) {
        double[][] profileMatrix = new double[4][motifs.get(0).length()];
        for (String dNa : motifs) {
            for (int i = 0; i < dNa.length(); i++) {
                profileMatrix[getNum(dNa.charAt(i))][i]++;
            }
        }

        if (applyLaplacian) {
            for (int i = 0; i < profileMatrix.length; i++) {
                for (int j = 0; j < profileMatrix[0].length; j++) {
                    profileMatrix[i][j]++;
                }
            }
        }

        return profileMatrix;
    }

    private int getNum(char c) {
        switch (c) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return -1;
        }
    }

    private char getProtien(int n) {
        switch (n) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        default:
            return '+';
        }
    }

    public static void main(String args[]) {
        In in = new In("src/GreedyMotifSearch");
        String nums[] = in.readLine().split(" ");
        int k = Integer.parseInt(nums[0]);
        int t = Integer.parseInt(nums[1]);

        List<String> dNaList = new ArrayList<String>();

        for (int i = 0; i < t; i++) {
            dNaList.add(in.readLine());
        }

        GreedyMotifSearch greedyMotifSearch = new GreedyMotifSearch(dNaList, k, t, true);
        List<String> bestMotifs = greedyMotifSearch.getBestMotifs();

        for (int i = 0; i < t; i++) {
            System.out.println(bestMotifs.get(i));
        }
    }
}
