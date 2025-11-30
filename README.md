  Ariel University Computer Science – Ex1
  Gal Attias  
  208743690


## Files in the Repository

- **Ex1.java**  
  Main class of the exercise. Implements all the required static methods on polynomials:​:contentReference[oaicite:3]{index=3}  
  - `f(double[] poly, double x)` – evaluates the polynomial at `x`.  
  - `root_rec(double[] p, double x1, double x2, double eps)` – recursively approximates a root in `[x1, x2]` using the bisection method.:contentReference[oaicite:4]{index=4}  
  - `PolynomFromPoints(double[] xx, double[] yy)` – builds a polynomial (up to degree 2) that interpolates 2 or 3 given points (returns `null` for invalid input).:contentReference[oaicite:5]{index=5}  
  - `equals(double[] p1, double[] p2)` – checks if two polynomials represent the same function up to `EPS`.  
  - `poly(double[] poly)` – converts a polynomial array to a human-readable string (e.g. `-1.2x^3 +3.1x^2 +2.0`).:contentReference[oaicite:6]{index=6}  
  - `getPolynomFromString(String p)` – parses a polynomial string back into a `double[]`, designed so that `getPolynomFromString(poly(p))` returns an array equal to `p`.:contentReference[oaicite:7]{index=7}  
  - `add(double[] p1, double[] p2)` – returns `p1 + p2`.  
  - `mul(double[] p1, double[] p2)` – returns the product polynomial `p1 * p2`.  
  - `derivative(double[] po)` – returns the derivative polynomial.  
  - `sameValue(double[] p1, double[] p2, double x1, double x2, double eps)` – finds an `x` such that `|p1(x) − p2(x)| < eps` using a recursive bisection-like method.:contentReference[oaicite:8]{index=8}  
  - `length(double[] p, double x1, double x2, int numberOfSegments)` – approximates the length of the curve `y = f(x)` between `x1` and `x2` using line segments.:contentReference[oaicite:9]{index=9}  
  - `area(double[] p1, double[] p2, double x1, double x2, int numberOfTrapezoid)` – approximates the area between two polynomials on `[x1, x2]` using a trapezoid-based Riemann-like method, including handling sign changes inside a segment.:contentReference[oaicite:10]{index=10}  

- **Ex1Test.java**  
  JUnit 5 test class for `Ex1`. It contains tests (given in the template and extended by me) that check:​:contentReference[oaicite:11]{index=11}  
  - `f` – evaluation of a polynomial at specific x values.  
  - `add` – properties such as commutativity and adding the negative polynomial.:contentReference[oaicite:12]{index=12}  
  - `mul` – multiplication with zero, commutativity, and consistency `p1(x) * p2(x) = (p1*p2)(x)`.:contentReference[oaicite:13]{index=13}  
  - `derivative` – repeated derivatives until reaching the zero polynomial.:contentReference[oaicite:14]{index=14}  
  - `poly` and `getPolynomFromString` – that parsing and printing are inverse operations.:contentReference[oaicite:15]{index=15}  
  - `equals` – equality and non-equality of different polynomial arrays.:contentReference[oaicite:16]{index=16}  
  - `sameValue` – symmetry of the function when swapping polynomials.:contentReference[oaicite:17]{index=17}  
  - `area` – symmetry (`area(p1,p2) = area(p2,p1)`) and correctness on simple functions like `0` and `x`.:contentReference[oaicite:18]{index=18}  

  **Additional tests I added:**:contentReference[oaicite:19]{index=19}  
  - `testAddZeroIdentity()` – checks that adding `ZERO` from the left or right does not change the polynomial.  
  - `testPolyAndFromString()` – checks that `poly` and `getPolynomFromString` behave as inverse operations on a polynomial with missing coefficients (zero terms in the middle).


## Ex1_GUI Output Image

<img width="503" height="505" alt="image" src="https://github.com/user-attachments/assets/4e805e85-0560-4355-a4c9-84f469b148c5" />

