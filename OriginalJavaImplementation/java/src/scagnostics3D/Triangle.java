package scagnostics3D;

public class Triangle {
	protected Point p1, p2, p3; 
	protected Edge e1, e2, e3; 
	// tetrahedra share this facet
	protected Tetrahedron T1 = null, T2 = null;	

	protected double area, perimeter;
	protected boolean onHull = false;
	protected boolean onShape = false;
	protected boolean isVisited = false;

	// normal to the plane determined  by this triangle
	protected Point normal ; 	
	// plane equation parameters. aX+bY+cZ+d=0, normal = (a,b,c)
	protected double a, b, c, d;  

	protected Triangle(Point n1, Point n2, Point n3){
		e1 = new Edge(n1, n2);
		e2 = new Edge(n2, n3);
		e3 = new Edge(n3, n1);
		update(e1, e2, e3);
	}

	protected Triangle( Edge a1, Edge a2, Edge a3 ) {
		update(a1, a2, a3);
	}

	public boolean Inside (Point x){
		return normal.Dot(x) > d;
	}

	protected void update(Edge e1, Edge e2, Edge e3) {
		e1.nextE = e2;
		e2.nextE = e3;
		e3.nextE = e1;

		p1 = e1.p1;
		p2 = e1.p2;
		if (e2.p1.isEqual(p2)){
			p3 = e2.p2;
		} else {
			p3 = e2.p1;
		}

		// compute the plane equation coefficients
		normal = p2.Subtract(p1).Cross(p3.Subtract(p1)).Normalize();

		a = normal.x;
		b = normal.y;
		c = normal.z;

		d = normal.Dot(p1);

		perimeter = e1.weight + e2.weight + e3.weight;
		double p = 0.5 * perimeter;
		area = Math.sqrt(p * (p - e1.weight) * (p - e2.weight) * (p - e3.weight));
		this.onHull = false;
		this.onShape = false;
	}

	protected boolean isEqual(Triangle t1){
		Point v1 = t1.p1;
		Point v2 = t1.p2;
		Point v3 = t1.p3;
		if (v1.isEqual(this.p1) && v2.isEqual(this.p2) && v3.isEqual(this.p3))
			return true;
		else
			return false;
	}

	protected boolean isEquivalent(Triangle t){
		Triangle t1 = new Triangle(t.p1, t.p3, t.p2);
		Triangle t2 = new Triangle(t.p2, t.p1, t.p3);
		Triangle t3 = new Triangle(t.p2, t.p3, t.p1);
		Triangle t4 = new Triangle(t.p3, t.p1, t.p2);
		Triangle t5 = new Triangle(t.p3, t.p2, t.p1);

		return ( this.isEqual(t) || this.isEqual(t1) || this.isEqual(t2) || this.isEqual(t3) 
				|| this.isEqual(t4) || this.isEqual(t5)); 
	}

}