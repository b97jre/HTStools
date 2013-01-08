package general;

public class SimpleLinearRegressionModel {
	/**
	* The predictor
	*/
	double[] x;
	
	/**
	* The response
	*/
	double[] y;
	
	/**
	* The mean of x
	*/
	double xBar;
	
	/**
	* The mean of y
	*/
	double yBar;
	
	/**
	* The slope
	*/
	double slope;
	
	/**
	* The intercept
	*/
	double intercept;
	
	/**
	* The std dev of x
	*/
	double sx;
	
	/**
	* The std dev of y
	*/
	double sy;
	
	
	/**
	* Constructor
	* Constructs a Simple Linear Regression Model (no checks are made on params)
	*
	* @param argX the predictor variable
	* @param argY the response variable
	*/ 
	public SimpleLinearRegressionModel(double[] argX, double[] argY) {
		x = argX;
		y = argY;
		compute();
	}
	
	private void compute() {
		int n = x.length-2;
		double sumy = 0.0,
		sumx = 0.0,
		sumx2 = 0.0,
		sumy2 = 0.0,
		sumxy = 0.0;
		
		for (int i = 0; i < n; i++) {
			sumx += x[i];
			sumx2 += x[i] * x[i];
			sumy += y[i];
			sumy2 += y[i] * y[i];
			sumxy += x[i] * y[i];
		}
		
		xBar = sumx / n;
		yBar = sumy / n;
		
		
		slope = (sumxy - sumx * yBar) / (sumx2 - sumx * xBar);
		intercept = yBar - slope * xBar;
		sx = Math.sqrt((sumx2 - sumx * xBar) / (n - 1));
		sy = Math.sqrt((sumy2 - sumy * yBar) / (n - 1));
	}
	
	/**
	* Get the slope for this model
	*
	* @return the slope
	*/
	public double getSlope() {
		return slope;
	}
	
	/**
	* Get the intercept for this model
	*
	* @return the intercept
	*/
	public double getIntercept() {
		return intercept;
	}
	
	/**
	* Get the R-squared for this model
	*
	* @return r-squared
	*/
	public double getRSquared() {
		double r = slope * sx / sy;
		return r * r;
	}
}



