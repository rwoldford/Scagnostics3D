package scagnostics3D;
import java.util.*;

class Point {
	protected double x, y, z;	// 3D coordinate X,Y,Z
	protected double[] v; 		// array of coordinates
	protected List<Edge> vE; 		// edges sharing this node 
	protected List<Triangle> vF;	// facets sharing this node
	protected List<Tetrahedron> vT; // tetrahedra sharing this node 
	protected boolean onMST;
	protected boolean isVisited = false;
	protected int mstDegree;
	protected int pointID;


	protected Point(double x, double y, double z, int pointID) {
		this.x = x;
		this.y = y;
		this.z = z;
		vE = new ArrayList<Edge>();
		vF = new ArrayList<Triangle>();
		vT = new ArrayList<Tetrahedron>();
		this.pointID = pointID;
		this.v = new double[3];
		v[0] = x;
		v[1] = y;
		v[2] = z;
	}

	protected Point(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
		vE = new ArrayList<Edge>();
		vF = new ArrayList<Triangle>();
		vT = new ArrayList<Tetrahedron>();
	}

	protected boolean isEqual(Point v){
		double zero = Double.MIN_VALUE;
		if ((Math.abs(v.x - this.x)<= zero) && (Math.abs(v.y - this.y)<=zero) 
				&& (Math.abs(v.z - this.z)<=zero))
			return true;
		else
			return false;

	}
	public Point Add(Point n){
		double ax = this.x + n.x;
		double ay = this.y + n.y;
		double az = this.z + n.z;

		return new Point(ax, ay, az);
	}

	public Point Subtract(Point n){
		double ax = this.x - n.x;
		double ay = this.y - n.y;
		double az = this.z - n.z;

		return new Point(ax, ay, az);
	}

	public Point Scale(double s){
		double ax = this.x * s;
		double ay = this.y * s;
		double az = this.z * s;
		return new Point(ax, ay, az);
	}

	public double Dot(Point n){
		return (this.x * n.x + this.y * n.y + this.z * n.z);
	}

	public Point Normalize(){
		return Scale( 1/Length() );
	}

	public double Length(){
		return Math.sqrt(Dot(this));
	}

	public Point Cross(Point n){
		return new Point((this.y * n.z - n.y * this.z),
				(this.z * n.x - n.z * this.x),
				(this.x * n.y - n.x * this.y) );
	}

	protected double distToPoint(double px, double py, double pz) {
		double dx = px - x;
		double dy = py - y;
		double dz = pz - z;
		return Math.sqrt(dx * dx + dy * dy + dz * dz);
	}


	protected Iterator<Edge> getNeighborIterator() {
		return vE.iterator();
	}

	protected Edge shortestEdge(boolean mst) {
		Edge emin = null;
		if (vE != null) {
			Iterator<Edge> it = vE.iterator();
			double wmin = Double.MAX_VALUE;
			while (it.hasNext()) {
				Edge e = (Edge) it.next();
				if (mst || !e.otherNode(this).onMST) {
					double wt = e.weight;
					if (wt < wmin) {
						wmin = wt;
						emin = e;
					}
				}
			}
		}
		return emin;
	}

	protected int getMSTChildren(double cutoff, double[] maxLength) {
		int count = 0;
		if (isVisited)
			return count;
		isVisited = true;
		Iterator<Edge> it = vE.iterator();
		while (it.hasNext()) {
			Edge e = (Edge) it.next();
			if (e.onMST) {
				if (e.weight < cutoff) {
					if (!e.otherNode(this).isVisited) {
						count += e.otherNode(this).getMSTChildren(cutoff, maxLength);
						double el = e.weight;
						if (el > maxLength[0])
							maxLength[0] = el;
					}
				}
			}
		}
		count++;
		return count;
	}
}