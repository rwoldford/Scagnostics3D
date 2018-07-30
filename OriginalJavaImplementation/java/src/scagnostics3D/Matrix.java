package scagnostics3D;

public class Matrix {
	protected int rows, cols;
	protected double[][] vMat;
	
public Matrix(double[][] pts){
	rows = pts.length + 1;
	cols = pts[0].length ;
	vMat = new double[rows][cols];

	if (cols > 3) {
		vMat[rows - 1] = new double[]{1, 1, 1, 1};
	} else {
		vMat[rows - 1] = new double[]{1, 1, 1};
	};
	for (int i = 0; i < rows - 1; i++ ){
		System.arraycopy(pts[i], 0, vMat[i], 0, cols);
	}
}

public Matrix(Point n1, Point n2, Point n3, Point n4){
	rows = 4 ;
	cols = 4 ;
    vMat = new double[4][4];
    vMat[0] = new double[]{n1.x, n2.x, n3.x, n4.x};
    vMat[1] = new double[]{n1.y, n2.y, n3.y, n4.y};
    vMat[2] = new double[]{n1.z, n2.z, n3.z, n4.z};
    vMat[3] = new double[]{1, 1, 1, 1};
}

public Matrix(double[] x, double[] y, double[] z) {
	rows = 3 ;
	cols = x.length ;
    vMat = new double[3][cols];
    vMat[0] = x;
    vMat[1] = y;
    vMat[2] = z;
}

public double Determinant(){
	return det(vMat);
}

private double det(double[][] mat) {

	double result = 0;

	if(mat.length == 1) {
		result = mat[0][0];
		return result;
	}

	if(mat.length == 2) {
		result = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
		return result;
	}

	for(int i = 0; i < mat[0].length; i++) {
		double temp[][] = new double[mat.length - 1][mat[0].length - 1];

		for(int j = 1; j < mat.length; j++) {
			for(int k = 0; k < mat[0].length; k++) {
				if(k < i) {
					temp[j - 1][k] = mat[j][k];
				} else if(k > i) {
					temp[j - 1][k - 1] = mat[j][k];
				} 
			}
		}
		result += mat[0][i] * Math.pow(-1, (double)i) * det(temp);
	}

	return result;
} 

public double[] Max() {
    
	double[] maximum = new double[rows]; // start with the first value
	double[] t = null;
	for (int r=0; r<rows; r++) {
		t = vMat[r];

		for (int i=0; i<t.length; i++) {
			if (t[i] > maximum[r]) {
				maximum[r] = t[i];
			}
		}
	}
    return maximum;
}

public double[] Min() {
    
	double[] minimum = new double[rows]; // start with the first value
	double[] t = null;
	for (int r = 0; r < rows; r++) {
		t = vMat[r];
		for (int i = 0; i < t.length; i++) {
			if (t[i] < minimum[r]) {
				minimum[r] = t[i];
			}
		}
	}
    return minimum;
}

}
