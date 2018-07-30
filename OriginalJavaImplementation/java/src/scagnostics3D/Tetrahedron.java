package scagnostics3D;

public class Tetrahedron {

	// viewed from n1, n2\n3\n4 in counter-clockwise order
	protected Point p1, p2, p3, p4;
	// facets of this tetrahedron
	protected Triangle F1, F2, F3, F4; 
    protected boolean onComplex = true;
    protected double volume; 

    protected Tetrahedron(Triangle t, Point p) {
        update(t, p);
    }

    protected Tetrahedron(Point p1, Point p2, Point p3, Point p4) {
    	Triangle t = new Triangle(p1, p2, p3);
        update(t, p4);
    }

    protected void update(Triangle t, Point p) {
        this.onComplex = true;
        this.p1 = t.p1;
        this.p2 = t.p2;
        this.p3 = t.p3;
        this.p4 = p;
        this.F1 = t;
        this.F2 = new Triangle(p1, p3, p4);
        this.F3 = new Triangle(p1, p4, p2);
        this.F4 = new Triangle(p3, p2, p4);
        Matrix Mat = new Matrix(p1, p2, p3, p4);
        double vol = Mat.Determinant()/6.0;
        if (vol > 0)
        {
        	this.p1 = t.p3;
        	this.p3 = t.p1;
        }
        this.volume = Math.abs(vol);
        this.onComplex = true;
    }
  
}