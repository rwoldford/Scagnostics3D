package scagnostics3D;

import java.util.Arrays;
import java.util.Comparator;

public class Sorts {

	private Sorts() {}

	public static void doubleArraySort(double[] x, int fromIndex, int toIndex) {
		if (fromIndex == toIndex) {
			fromIndex = 0;
			toIndex = x.length;
		}
		Arrays.sort(x, fromIndex, toIndex);
	}

	@SuppressWarnings("unchecked")
	public static int[] indexedDoubleArraySort(final double[] x, int fromIndex, int toIndex) {
		if (fromIndex == toIndex) {
			fromIndex = 0;
			toIndex = x.length;
		}
		Integer[] sortOrder = new Integer[toIndex - fromIndex];
		for (int i = 0; i < sortOrder.length; i++) {
			sortOrder[i] = new Integer(i);
		}
		Arrays.sort(sortOrder, fromIndex, toIndex, new Comparator() {
			public int compare(Object object1, Object object2) {
				int firstIndex = ((Integer) object1).intValue();
				int secondIndex = ((Integer) object2).intValue();
				return Double.compare(x[firstIndex], x[secondIndex]);
			}
		});
		int[] result = new int[sortOrder.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = sortOrder[i].intValue();
		}
		return result;
	}

	public static double[] rank(double[] a) {

		int k, k1, k2, kind, kms, l, lind, n;
		double ak, am, freq;
		boolean insert;

		n = a.length;

		double[] ranks = new double[n];

		int[] index = indexedDoubleArraySort(a, 0, n);

		lind = index[0];
		am = a[lind];
		k1 = 0;
		k2 = 0;
		ak = 1.0;
		/* kms allows for missing data */
		kms = 1;
		for (k = 1; k < n; k++) {
			kind = index[k];
			insert = true;
			if (!Double.isNaN(am)) {
				freq = 1.0;
				/*
                if (wt != null)
                    freq = Math.floor(wt[kind]);
				 */
				kms += (int) freq;
				if (freq > 1.0) {
					ak += 0.5 * (freq - 1.0);
					k1 = k;
					k2 = k;
				} else if (a[kind] == am) {
					k2 = k;
					ak += 0.5;
					if (k < n - 1)
						insert = false;
				}
				if (insert) {
					for (l = k1; l <= k2; l++) {
						lind = index[l];
						ranks[lind] = ak;
					}
					if (k2 != n - 1 && k == n - 1)
						ranks[kind] = kms;
				}
			}
			if (insert) {
				k1 = k;
				k2 = k;
				ak = kms;
				am = a[kind];
			}
		}
		return ranks;
	}
}
