package scagnostics3D;

import java.io.File;
import java.util.*;

public class Main {
	static String osname = System.getProperty("os.name");

    public static void main(String[] argv) {

    	String fname = "/tmp/scagnostics3D/iris.txt";

    	if (osname.equalsIgnoreCase("Windows XP")){
    		fname = "C:\\Test\\iris.txt";
    	}
    	
    	double[][] points = getData(new File(fname));
    	double[][] scagnostics = computeScagnostics(points);
        // test
        String[] sl = Scagnostics3D.getScagnosticsLabels();
        for (int l=0;l<scagnostics.length;l++){
        	for (int kk=0; kk<9; kk++)
        		System.out.println(sl[kk] + ": " + scagnostics[l][kk]);
        	System.out.println("---------------------- ");
        }
        boolean[] outliers = Scagnostics3D.computeScagnosticsOutliers(scagnostics);
        boolean[] exemplars = Scagnostics3D.computeScagnosticsExemplars(scagnostics);
        computeTests(scagnostics, outliers, exemplars);
    }

    private static double[][] getData(File fname) {
        java.io.BufferedReader fin;
        try {
            fin = new java.io.BufferedReader(new java.io.FileReader(fname));
        } catch (java.io.FileNotFoundException fe) {
            javax.swing.JOptionPane.showMessageDialog(null, "File not found!", "Alert",
                    javax.swing.JOptionPane.ERROR_MESSAGE);
            return null;
        }
        try {
            String record;
            record = fin.readLine();
            //Get the scagnosticsLabels
            record = replaceSeparatorsWithBlanks(record);
            StringTokenizer st = new StringTokenizer(record, " ");
            int col = 0;
            int numVars = st.countTokens();
            String[] variableLabels = new String[numVars];

            while (st.hasMoreTokens()) {
                variableLabels[col] = st.nextToken();
                col++;
            }
            //Count the number of rows
            int numRows = 0;
            record = fin.readLine();
            while (record != null) {
                record = fin.readLine();
                numRows++;
            }
            fin.close();
            System.out.println("Number of rows, cols " + numRows + " " + numVars);

            //Read in the data
            fin = new java.io.BufferedReader(new java.io.FileReader(fname));
            double[][] data = new double[numVars][numRows];
            record = fin.readLine();    //ignore line with scagnosticsLabels
            record = fin.readLine();
            int j = 0;
            while (record != null) {
                record = replaceSeparatorsWithBlanks(record);
                st = new StringTokenizer(record, " ");
                for (int i = 0; i < numVars; i++) {
                    try {
                        String tmp = st.nextToken();
                        data[i][j] = Double.parseDouble(tmp);
                    } catch (Exception ie) {
                        javax.swing.JOptionPane.showMessageDialog(null, 
                        		"Error reading from the file", "Alert",
                                javax.swing.JOptionPane.ERROR_MESSAGE);
                        return null;
                    }
                }
                record = fin.readLine();
                j++;
            }
            fin.close();

            return data;
        } catch (java.io.IOException ie) {
            javax.swing.JOptionPane.showMessageDialog(null, 
            		"Error reading from the file", "Alert",
                    javax.swing.JOptionPane.ERROR_MESSAGE);
            return null;
        }
    }

    private static String replaceSeparatorsWithBlanks(String record) {
        record = replaceAll(record, ",,", ",~,");
        record = replaceAll(record, "\t\t", "\t~\t");
        record = replaceAll(record, ",", " ");
        record = replaceAll(record, "\t", " ");
        return record;
    }

    private static String replaceAll(String source, String toReplace, String replacement) {
        int idx = source.lastIndexOf(toReplace);
        if (idx != -1) {
            StringBuffer sb = new StringBuffer(source);
            sb.replace(idx, idx + toReplace.length(), replacement);
            while ((idx = source.lastIndexOf(toReplace, idx - 1)) != -1) {
                sb.replace(idx, idx + toReplace.length(), replacement);
            }
            source = sb.toString();
        }
        return source;
    }

    private static double[][] computeScagnostics(double[][] points ) {
        normalizePoints(points);
        int nDim = points.length;
        int numCells = nDim * (nDim - 1)* (nDim - 2) / 6;
        double[][] scagnostics = new double[numCells][Scagnostics3D.getNumScagnostics()];
        int l = 0;
        for (int i = 2; i < nDim; i++) {
            for (int j = 1; j < i; j++) {
            	for (int k = 0; k < j; k++) {
                Scagnostics3D s = new Scagnostics3D(points[k], points[j], points[i]);
                scagnostics[l] = s.compute();
                l++;
            	}
            }
        }
        return scagnostics;
    }

    private static void normalizePoints(double[][] points) {
        double[] min = new double[points.length];
        double[] max = new double[points.length];
        for (int i = 0; i < points.length; i++) {
            min[i] = Double.MAX_VALUE;
            max[i] = Double.MIN_VALUE;
        }
        for (int i = 0; i < points.length; i++) {
            for (int j = 0; j < points[0].length; j++) {
                if (min[i] > points[i][j]) min[i] = points[i][j];
                if (max[i] < points[i][j]) max[i] = points[i][j];
            }
        }
        for (int i = 0; i < points.length; i++) {
            for (int j = 0; j < points[0].length; j++) {
                points[i][j] = (points[i][j] - min[i]) / (max[i] - min[i]);
            }
        }
    }

    private static void computeTests(double[][] scagnostics, boolean[] outliers, boolean[] exemplars) {
        double[][] test = getData(new File("scagnostics.txt"));

        for (int i = 0; i < test.length; i++) {
            for (int j = 0; j < test[0].length; j++) {
                if (scagnostics[j][i] != test[i][j]) System.out.println("error " + i + " " + j);
            }
        }

        for (int i = 0; i < outliers.length; i++) {
            if (i == 1 || i == 2 || i == 34 || i == 90) {
                if (!outliers[i]) System.out.println("error " + i);
            } else {
                if (outliers[i]) System.out.println("error " + i);
            }
        }

        for (int i = 0; i < exemplars.length; i++) {
            if (i == 13 || i == 95 || i == 118) {
                if (!exemplars[i]) System.out.println("error " + i);
            } else {
                if (exemplars[i]) System.out.println("error " + i);
            }
        }
    }
}