import java.util.ArrayList;
import java.util.List;
import java.util.Stack;
import java.util.TreeSet;

import edu.princeton.cs.algs4.In;

public class Bioinformatics {

    public static int computeHamming(String first, String second) {
        if (first == null || second == null) {
            throw new IllegalArgumentException("One or both the strings are null");
        }

        if (first.length() != second.length()) {
            throw new IllegalArgumentException("String lengths are not same");
        }

        int hammingDistance = 0;

        for (int i = 0; i < first.length(); i++) {
            if (first.charAt(i) != second.charAt(i)) {
                hammingDistance = hammingDistance + 1;
            }
        }

        return hammingDistance;
    }
    
    public static int count(String text, String pattern, int k) {
        int count = 0;
        int m = pattern.length();
        for (int i = 0; i <= text.length() - m ; i++) {
            if (computeHamming(text.substring(i, i + m), pattern) <= k) {
                count = count + 1;
            }
        }
        
        return count;
    }

    public static int minSkew(String dNa) {
        int gCount = 0;
        int cCount = 0;

        int minSkew = Integer.MAX_VALUE;

        for (int i = 0; i < dNa.length(); i++) {
            if (dNa.charAt(i) == 'G') {
                gCount = gCount + 1;
            } else if (dNa.charAt(i) == 'C') {
                cCount = cCount + 1;
            }

            int currSkew = gCount - cCount;
            if (currSkew < minSkew) {
                minSkew = currSkew;
            }
        }

        return minSkew;
    }

    public static int maxSkew(String dNa) {
        int gCount = 0;
        int cCount = 0;

        int maxSkew = Integer.MIN_VALUE;

        for (int i = 0; i < dNa.length(); i++) {
            if (dNa.charAt(i) == 'G') {
                gCount = gCount + 1;
            } else if (dNa.charAt(i) == 'C') {
                cCount = cCount + 1;
            }

            int currSkew = gCount - cCount;
            if (currSkew > maxSkew) {
                maxSkew = currSkew;
            }
        }

        return maxSkew;
    }

    public static ArrayList<Integer> minSkewPos(String dNa) {
        int gCount = 0;
        int cCount = 0;

        ArrayList<Integer> minSkewPositions = new ArrayList<Integer>();
        int skewValues[] = new int[dNa.length()];
        int minSkew = Integer.MAX_VALUE;

        for (int i = 0; i < dNa.length(); i++) {
            if (dNa.charAt(i) == 'G') {
                gCount = gCount + 1;
            } else if (dNa.charAt(i) == 'C') {
                cCount = cCount + 1;
            }

            int currSkew = gCount - cCount;
            skewValues[i] = currSkew;
            if (currSkew < minSkew) {
                minSkew = currSkew;
            }
        }

        for (int i = 0; i < dNa.length(); i++) {
            if (skewValues[i] == minSkew) {
                minSkewPositions.add(i + 1);
            }
        }

        return minSkewPositions;
    }
    
    public static ArrayList<Integer> maxSkewPos(String dNa) {
        int gCount = 0;
        int cCount = 0;

        ArrayList<Integer> maxSkewPositions = new ArrayList<Integer>();
        int skewValues[] = new int[dNa.length()];
        int maxSkew = Integer.MIN_VALUE;

        for (int i = 0; i < dNa.length(); i++) {
            if (dNa.charAt(i) == 'G') {
                gCount = gCount + 1;
            } else if (dNa.charAt(i) == 'C') {
                cCount = cCount + 1;
            }

            int currSkew = gCount - cCount;
            skewValues[i] = currSkew;
            if (currSkew > maxSkew) {
                maxSkew = currSkew;
            }
        }

        for (int i = 0; i < dNa.length(); i++) {
            if (skewValues[i] == maxSkew) {
                maxSkewPositions.add(i + 1);
            }
        }

        return maxSkewPositions;
    }

    public static ArrayList<Integer> approximateMatching(String text, String pattern, int k) {
        ArrayList<Integer> approximateMatchingPositions = new ArrayList<Integer>();

        int m = pattern.length();
        for (int i = 0; i <= text.length() - m; i++) {
            if (computeHamming(text.substring(i, i + m), pattern) <= k) {
                approximateMatchingPositions.add(i);
            }
        }

        return approximateMatchingPositions;
    }

    public static ArrayList<String> generateAllPermutations(int k) {
        ArrayList<String> permutations = new ArrayList<String>();
        generateAllPermutations(new Stack<Character>(), k, permutations);
        return permutations;
    }

    private static void generateAllPermutations(Stack<Character> stack, int k, ArrayList<String> permutations) {
        if (k == 0) {
            StringBuilder str = new StringBuilder();
            for (Character ch : stack) {
                str.append(ch);
            }
            permutations.add(str.toString());

            return;
        }

        stack.push('A');
        generateAllPermutations(stack, k - 1, permutations);
        stack.pop();

        stack.push('C');
        generateAllPermutations(stack, k - 1, permutations);
        stack.pop();

        stack.push('G');
        generateAllPermutations(stack, k - 1, permutations);
        stack.pop();

        stack.push('T');
        generateAllPermutations(stack, k - 1, permutations);
        stack.pop();
    }

    public static ArrayList<String> frequentKMerWithMismatches(String text, int k, int d) {
        ArrayList<String> frequentKMersWithMismatches = new ArrayList<String>();
        int[] freq = new int[(int) Math.pow(4, k)];

        int maxVal = Integer.MIN_VALUE;

        int index = 0;
        ArrayList<String> allPermutations = generateAllPermutations(k);
        for (String pattern : allPermutations) {
            int len = approximateMatching(text, pattern, d).size();
            freq[index] = len;

            if (len > maxVal) {
                maxVal = len;
            }

            index = index + 1;
        }

        for (int i = 0; i < freq.length; i++) {
            if (freq[i] == maxVal) {
                frequentKMersWithMismatches.add(allPermutations.get(i));
            }
        }

        return frequentKMersWithMismatches;
    }

    public static ArrayList<String> frequentKMerWithMismatchesReverse(String text, int k, int d) {
        ArrayList<String> frequentKMersWithMismatches = new ArrayList<String>();
        int[] freq = new int[(int) Math.pow(4, k)];

        int maxVal = Integer.MIN_VALUE;

        int index = 0;
        ArrayList<String> allPermutations = generateAllPermutations(k);
        for (String pattern : allPermutations) {
            int len = approximateMatching(text, pattern, d).size()
                    + approximateMatching(text, getReverseComplement(pattern), d).size();
            freq[index] = len;

            if (len > maxVal) {
                maxVal = len;
            }

            index = index + 1;
        }

        for (int i = 0; i < freq.length; i++) {
            if (freq[i] == maxVal) {
                frequentKMersWithMismatches.add(allPermutations.get(i));
            }
        }

        return frequentKMersWithMismatches;
    }

    public static String getReverseComplement(String protien) {
        return getComplement(new StringBuilder(protien).reverse().toString());
    }

    public static String getComplement(String protien) {
        StringBuilder complement = new StringBuilder();
        for (int i = 0; i < protien.length(); i++) {
            complement.append(getComplement(protien.charAt(i)));
        }

        return complement.toString();
    }

    private static char getComplement(char ch) {
        switch (ch) {
        case 'A':
            return 'T';
        case 'T':
            return 'A';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        default:
            return '+';
        }
    }

    public static ArrayList<String> motifEnumeration(List<String> dNa, int k, int d) {
        TreeSet<String> patterns = new TreeSet<String>();

        String processedDna = preprocessDna(dNa);
        for (int i = 0; i <= processedDna.length() - k; i++) {
            String currKMer = processedDna.substring(i, i + k);

            if (currKMer.equals("ATA"))
                System.out.println("currKMer : " + currKMer);
            for (String mismatchedPattern : neighbors(currKMer, d)) {
                if (patternAppearsInAllStringWithMismatches(dNa, mismatchedPattern, d)) {
                    patterns.add(mismatchedPattern);
                }
            }
        }

        return new ArrayList<String>(patterns);
    }

    private static boolean patternAppearsInAllStringWithMismatches(List<String> dNa, String mismatchedPattern, int d) {
        for (String pattern : dNa) {
            if (approximateMatching(pattern, mismatchedPattern, d).isEmpty()) {
                return false;
            }
        }

        return true;
    }

    public static long patternToNumber(String pattern) {
        int patternLen = pattern.length();
        long num = 0;
        long multiplier = 1;

        for (int i = 0; i < patternLen; i++) {
            num = num + ((long) getNum(pattern.charAt(patternLen - i - 1))) * multiplier;
            multiplier = multiplier << 2;
        }

        return num;
    }

    private static int getNum(char c) {
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

    public static String numberToPattern(int n, int k) {
        StringBuilder pattern = new StringBuilder();
        while (n != 0) {
            pattern.append(getProtien(n % 4));
            n = n / 4;
        }

        for (int i = 0; i < k - pattern.length(); i++) {
            pattern.append('A');
        }

        return pattern.reverse().toString();
    }

    private static char getProtien(int n) {
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

    public static ArrayList<String> neighbors(String currKMer, int d) {
        ArrayList<String> neighbors = new ArrayList<String>();
        for (String pattern : generateAllPermutations(currKMer.length())) {
            if (computeHamming(pattern, currKMer) <= d) {
                neighbors.add(pattern);
            }
        }

        return neighbors;
    }

    private static String preprocessDna(List<String> dNa) {
        StringBuilder combinedDna = new StringBuilder();
        for (String protien : dNa) {
            combinedDna.append(protien);
        }

        return combinedDna.toString();
    }

    public static String medianString(int k, List<String> dNaList) {
        int minScore = Integer.MAX_VALUE;
        String bestPattern = "";

        for (String pattern : generateAllPermutations(k)) {
            int currScore = getMinScore(dNaList, pattern);
            if (currScore < minScore) {
                minScore = currScore;
                bestPattern = pattern;
            }
        }

        return bestPattern;
    }

    public static int getMinScore(List<String> dNaList, String pattern) {
        int k = pattern.length();
        int score = 0;
        for (String dNa : dNaList) {
            int minHamm = Integer.MAX_VALUE;
            for (int i = 0; i <= dNa.length() - k; i++) {
                int currHamm = computeHamming(pattern, dNa.substring(i, i + k));
                minHamm = Math.min(minHamm, currHamm);
            }

            score = score + minHamm;
        }

        return score;
    }

    public static String profileMostProbable(String text, int k, double profile[][]) {
        double maxProb = -1.0;;
        String maxProbPattern = "";
        for (int i = 0; i <= text.length() - k; i++) {
            double currProb = 1.00;
            String currPatt = text.substring(i, i + k);
            for (int j = 0; j < currPatt.length(); j++) {
                currProb *= profile[(int) patternToNumber("" + currPatt.charAt(j))][j];
            }

            if (currProb > maxProb) {
                maxProb = currProb;
                maxProbPattern = currPatt;
            }
        }

        return maxProbPattern;
    }

    public static int[] computeFrequencies(String text, int k) {
        long num[] = new long[text.length()];
        int freq[] = new int[(int) Math.pow(4, k)];

        num[k - 1] = patternToNumber(text.substring(0, k));
        int multiplier = (int) Math.pow(4, k - 1);

        for (int i = k; i < text.length(); i++) {
            num[i] = (num[i - 1] - multiplier * getNum(text.charAt(i - k))) * 4 + getNum(text.charAt(i));
        }

        for (int i = k - 1; i < text.length(); i++) {
            freq[(int) num[i]]++;
        }

        return freq;
    }

    public static void main(String args[]) {
        In in = new In("src/minimum_skew_dataset");
        String text[] = in.readLine().split(" ");
        int k = Integer.parseInt(text[0]);
        // int d = Integer.parseInt(text[1]);

        List<String> list = new ArrayList<String>();
        while (!in.isEmpty()) {
            list.add(in.readLine());
        }

        System.out.println(medianString(k, list));
        System.out.println(count("CGTGACAGTGTATGGGCATCTTT", "TGT", 1));
    }
}
