

/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 *
 * This class represents a set of static methods on polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main class you should implement (see "add your code below").
 */
public class Ex1 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynomial function is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};

	/**
	 * Computes the f(x) value of the polynomial function at x.
	 * @param poly - polynomial function
	 * @param x value
	 * @return f(x) - the polynomial function value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for (int i = 0; i < poly.length; i++) {
			double c = Math.pow(x, i);
			ans += c * poly[i];
		}
		return ans;
	}

	/**
	 * Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps,
	 * assuming p(x1)*p(x2) <= 0.
	 * This function is implemented recursively using the bisection method.
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p, x1);
		double x12 = (x1 + x2) / 2.0;
		double f12 = f(p, x12);
		if (Math.abs(f12) < eps) {
			return x12;
		}
		if (f12 * f1 <= 0) {
			return root_rec(p, x1, x12, eps);
		} else {
			return root_rec(p, x12, x2, eps);
		}
	}

	/**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx x coordinates
	 * @param yy y coordinates
	 * @return an array of doubles representing the coefficients of the polynom.
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double[] ans = null;
		if (xx == null || yy == null) {
			return null;
		}
		int lx = xx.length;
		int ly = yy.length;
		if (lx != ly || lx <= 1 || lx >= 4) {
			return null;
		}

		// 2 points -> line: a0 + a1 x
		if (lx == 2) {
			double x0 = xx[0], y0 = yy[0];
			double x1 = xx[1], y1 = yy[1];
			double a1 = (y1 - y0) / (x1 - x0);
			double a0 = y0 - a1 * x0;
			ans = new double[2];
			ans[0] = a0;
			ans[1] = a1;
		}
		// 3 points -> degree <= 2, using Lagrange interpolation
		else { // lx == 3
			ans = new double[3]; // all zeros by default
			for (int i = 0; i < 3; i++) {
				double xi = xx[i];
				double yi = yy[i];

				// start with polynomial 1
				double[] li = {1.0};

				for (int j = 0; j < 3; j++) {
					if (j == i) continue;
					double xj = xx[j];
					double denom = xi - xj;

					// multiply current li by (x - xj)/denom = (-xj/denom) + (1/denom)*x
					double a0 = -xj / denom;
					double a1 = 1.0 / denom;
					double[] tmp = new double[li.length + 1];
					for (int k = 0; k < li.length; k++) {
						tmp[k] += li[k] * a0;
						tmp[k + 1] += li[k] * a1;
					}
					li = tmp;
				}

				// add yi * li to answer
				for (int k = 0; k < li.length; k++) {
					ans[k] += yi * li[k];
				}
			}
		}
		return ans;
	}

	/**
	 * Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
	 * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
	 * In practice here we normalize trailing exact zeros and compare coefficients up to EPS.
	 * @param p1 first polynomial function
	 * @param p2 second polynomial function
	 * @return true iff p1 represents the same polynomial function as p2.
	 */
	public static boolean equals(double[] p1, double[] p2) {
		if (p1 == p2) {
			return true;
		}
		if (p1 == null || p2 == null) {
			return false;
		}
		double[] a = trimZerosExact(p1);
		double[] b = trimZerosExact(p2);

		if (a.length != b.length) {
			return false;
		}
		for (int i = 0; i < a.length; i++) {
			if (Math.abs(a[i] - b[i]) > EPS) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Computes a String representing the polynomial function.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynomial function represented as an array of doubles
	 * @return String representing the polynomial function.
	 */
	public static String poly(double[] poly) {
		if (poly == null || poly.length == 0) {
			return "0";
		}
		double[] p = trimZerosExact(poly);
		if (p.length == 1 && p[0] == 0.0) {
			return "0";
		}

		StringBuilder sb = new StringBuilder();
		int maxDeg = p.length - 1;
		boolean first = true;

		for (int i = maxDeg; i >= 0; i--) {
			double c = p[i];
			if (Math.abs(c) <= EPS) {
				continue;
			}
			double absC = Math.abs(c);

			if (first) {
				if (c < 0) {
					sb.append("-");
				}
				first = false;
			} else {
				if (c < 0) {
					sb.append(" -");
				} else {
					sb.append(" +");
				}
			}

			sb.append(absC);
			if (i >= 1) {
				sb.append("x");
				if (i > 1) {
					sb.append("^").append(i);
				}
			}
		}

		String ans = sb.toString();
		if (ans.length() == 0) {
			ans = "0";
		}
		return ans;
	}

	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps.
	 * This function computes an x value (x1<=x<=x2) for which |p1(x) -p2(x)| < eps,
	 * assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * The implementation is similar to root_rec, applied to the difference p1-p2.
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
		double f1 = f(p1, x1) - f(p2, x1);
		double mid = (x1 + x2) / 2.0;
		double fm = f(p1, mid) - f(p2, mid);

		if (Math.abs(fm) < eps || Math.abs(x2 - x1) < eps) {
			return mid;
		}
		if (f1 * fm <= 0) {
			return sameValue(p1, p2, x1, mid, eps);
		} else {
			return sameValue(p1, p2, mid, x2, eps);
		}
	}

	/**
	 * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
	 * This function computes an approximation of the length of the function between f(x1) and f(x2)
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2.
	 * This function is implemented iteratively (non recursive).
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
		if (p == null || numberOfSegments <= 0) {
			return 0;
		}
		double ans = 0;
		double dx = (x2 - x1) / numberOfSegments;

		double xPrev = x1;
		double yPrev = f(p, xPrev);

		for (int i = 1; i <= numberOfSegments; i++) {
			double x = x1 + i * dx;
			double y = f(p, x);
			double dxSeg = x - xPrev;
			double dySeg = y - yPrev;
			ans += Math.sqrt(dxSeg * dxSeg + dySeg * dySeg);
			xPrev = x;
			yPrev = y;
		}
		return ans;
	}

	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer n which represents the number of
	 * trapezoids between the functions (number of samples on each polynom).
	 * This function computes an approximation of the area between the polynomial functions within the x-range.
	 * The area is computed using a Riemann-like integral.
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of trapezoids between x1 and x2.
	 * @return the approximated area between the two polynomial functions within the [x1,x2] range.
	 */
	public static double area(double[] p1, double[] p2, double x1, double x2, int numberOfTrapezoid) {
		if (p1 == null || p2 == null || numberOfTrapezoid <= 0) {
			return 0;
		}
		double ans = 0;
		double dx = (x2 - x1) / numberOfTrapezoid;

		for (int i = 0; i < numberOfTrapezoid; i++) {
			double xa = x1 + i * dx;
			double xb = xa + dx;
			double ya = f(p1, xa) - f(p2, xa);
			double yb = f(p1, xb) - f(p2, xb);

			// if same sign: simple trapezoid on |difference|
			if (ya * yb >= 0) {
				ans += 0.5 * (Math.abs(ya) + Math.abs(yb)) * dx;
			} else {
				// sign change inside [xa, xb] -> split at intersection point
				double xm = sameValue(p1, p2, xa, xb, EPS);
				double ym = 0.0; // at the intersection the difference is 0

				ans += 0.5 * (Math.abs(ya) + Math.abs(ym)) * (xm - xa);
				ans += 0.5 * (Math.abs(ym) + Math.abs(yb)) * (xb - xm);
			}
		}
		return ans;
	}

	/**
	 * This function computes the array representation of a polynomial function from a String
	 * representation. Note: given a polynomial function represented as a double array,
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * @param p - a String representing polynomial function.
	 * @return polynomial as array of doubles.
	 */
	public static double[] getPolynomFromString(String p) {
		if (p == null) {
			return ZERO;
		}
		String s = p.replace(" ", "").toLowerCase();
		if (s.length() == 0) {
			return ZERO;
		}
		if (s.charAt(0) != '+' && s.charAt(0) != '-') {
			s = "+" + s;
		}

		// first pass: find max degree
		int i = 0;
		int maxDeg = 0;
		while (i < s.length()) {
			// sign
			i++; // skip '+' or '-'
			int j = i;
			while (j < s.length() && s.charAt(j) != '+' && s.charAt(j) != '-') {
				j++;
			}
			String term = s.substring(i, j);
			int xIndex = term.indexOf('x');
			int deg;
			if (xIndex == -1) {
				deg = 0;
			} else {
				if (xIndex == term.length() - 1) {
					deg = 1;
				} else {
					int powIndex = term.indexOf('^', xIndex);
					if (powIndex == -1) {
						deg = 1;
					} else {
						String powStr = term.substring(powIndex + 1);
						deg = Integer.parseInt(powStr);
					}
				}
			}
			if (deg > maxDeg) {
				maxDeg = deg;
			}
			i = j;
		}

		double[] ans = new double[maxDeg + 1];

		// second pass: fill coefficients
		i = 0;
		while (i < s.length()) {
			char sign = s.charAt(i);
			int signVal = (sign == '-') ? -1 : 1;
			i++;
			int j = i;
			while (j < s.length() && s.charAt(j) != '+' && s.charAt(j) != '-') {
				j++;
			}
			String term = s.substring(i, j);
			int xIndex = term.indexOf('x');
			double coef;
			int deg;

			if (xIndex == -1) { // constant term
				coef = Double.parseDouble(term);
				deg = 0;
			} else {
				String coefPart = term.substring(0, xIndex);
				if (coefPart.length() == 0) {
					coef = 1.0;
				} else {
					coef = Double.parseDouble(coefPart);
				}
				if (xIndex == term.length() - 1) {
					deg = 1;
				} else {
					int powIndex = term.indexOf('^', xIndex);
					if (powIndex == -1) {
						deg = 1;
					} else {
						String powStr = term.substring(powIndex + 1);
						deg = Integer.parseInt(powStr);
					}
				}
			}

			ans[deg] += signVal * coef;
			i = j;
		}

		ans = trimZerosExact(ans);
		return ans;
	}

	/**
	 * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
	 * @param p1 first polynomial
	 * @param p2 second polynomial
	 * @return p1 + p2
	 */
	public static double[] add(double[] p1, double[] p2) {
		if (p1 == null && p2 == null) {
			return ZERO;
		}
		if (p1 == null) {
			return trimZerosExact(p2.clone());
		}
		if (p2 == null) {
			return trimZerosExact(p1.clone());
		}

		int n = Math.max(p1.length, p2.length);
		double[] ans = new double[n];
		for (int i = 0; i < n; i++) {
			double a = (i < p1.length) ? p1[i] : 0.0;
			double b = (i < p2.length) ? p2[i] : 0.0;
			ans[i] = a + b;
		}
		ans = trimZerosExact(ans);
		return ans;
	}

	/**
	 * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
	 * @param p1 first polynomial
	 * @param p2 second polynomial
	 * @return p1 * p2
	 */
	public static double[] mul(double[] p1, double[] p2) {
		if (p1 == null || p2 == null) {
			return ZERO;
		}
		double[] a = trimZerosExact(p1);
		double[] b = trimZerosExact(p2);
		if (a.length == 1 && a[0] == 0.0) {
			return ZERO;
		}
		if (b.length == 1 && b[0] == 0.0) {
			return ZERO;
		}

		double[] ans = new double[a.length + b.length - 1];
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < b.length; j++) {
				ans[i + j] += a[i] * b[j];
			}
		}
		ans = trimZerosExact(ans);
		return ans;
	}

	/**
	 * This function computes the derivative of the p0 polynomial function.
	 * @param po polynomial
	 * @return derivative polynomial
	 */
	public static double[] derivative(double[] po) {
		if (po == null || po.length == 0) {
			return ZERO;
		}
		double[] p = trimZerosExact(po);
		if (p.length <= 1) {
			return ZERO;
		}
		double[] ans = new double[p.length - 1];
		for (int i = 1; i < p.length; i++) {
			ans[i - 1] = i * p[i];
		}
		ans = trimZerosExact(ans);
		return ans;
	}

	/**
	 * Helper: removes exact trailing zeros from a polynomial representation.
	 * Keeps at least one coefficient (so {0,0} -> {0}).
	 */
	private static double[] trimZerosExact(double[] p) {
		if (p == null || p.length == 0) {
			return ZERO;
		}
		int last = p.length - 1;
		while (last > 0 && p[last] == 0.0) {
			last--;
		}
		double[] ans = new double[last + 1];
		for (int i = 0; i <= last; i++) {
			ans[i] = p[i];
		}
		return ans;
	}

}
