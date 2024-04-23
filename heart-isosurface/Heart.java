
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;
import java.lang.Math;

class Volume
{
    // Data structure to store sampled heart equation
	int data[][][];
	int imgResolution = 512;
	int gridResolution = 256;
    float zoom = 1.f;
	float samplingDistance=0.5f;

	int octreeBlockSize = 8;
	public int octree[][][] = null; // 8x8x8 voxel blocks 


	// Getters and Setters
	// -------------------
	public int getImageResolution() 
	{
		return imgResolution;
	}
	
	public void setImageResolution(int imgReso) 
	{
		imgResolution = imgReso;
	}
	
	public int getGridResolution() 
	{
		return gridResolution;
	}
	
	public void setGridResolution(int gridReso) 
	{
		gridResolution = gridReso;
	}

	public float getZoom() {
		return zoom;
	}

	public void setZoom() {
		this.zoom = (float) imgResolution / (float) gridResolution;
	}

	public float getSamplingDistance() {
		return samplingDistance;
	}

	public void setSamplingDistance(float samplingDistance) {
		this.samplingDistance = samplingDistance;
	}
	// -------------------

    // Convert value from one range to another
    double range(double value, double low, double high, double rlow, double rhigh) {
        return rlow + ((rhigh - rlow) * ((value - low) / (high - low)));
    }

	// Sample the heart equation and store it in a 3D array
    void makeHeart() {
        data = new int[gridResolution][gridResolution][gridResolution];
		int [][][] newvol=new int[gridResolution][gridResolution][gridResolution];
        for (int x = 0; x < gridResolution; x++) {
			double xScaled = range(x, 0, gridResolution-1, -2, 2); //= -2 + x * 4.0 / (gridResolution - 1); //scale z to -2 to 2

            for (int y = 0; y < gridResolution; y++) {
				double yScaled = range(y, 0, gridResolution-1, -2, 2); //= -2 + y * 4.0 / (gridResolution - 1); //scale y to -2 to 2

                for (int z = 0; z < gridResolution; z++) {
					double zScaled = range(z, 0, gridResolution-1, -2, 2);  //= -2 + z * 4.0 / (gridResolution - 1); //scale x to -2 to 2

                    double ans = (- (Math.pow((Math.pow(xScaled, 2) + 
													2 * Math.pow(yScaled, 2) + 
													Math.pow(zScaled, 2) - 1), 3) -
												Math.pow(xScaled, 2) * Math.pow(zScaled, 3) - 
												0.1 * Math.pow(yScaled, 2) * Math.pow(zScaled, 3)));
					//scale ans to 0 to 255: [0, 0.00001] -> [0, 255], negative values get scaled too but remain negative
					newvol[x][y][z] = (int) range(ans, 0, 0.00001, 0, gridResolution - 1); // ((gridResolution - 1.0) * ans/0.00001); 
                }
            }
        }
		// to rotate the volume
		for (int k = 0; k < gridResolution; k++) 
			for (int j = 0; j < gridResolution; j++) 
				for (int i = 0; i < gridResolution; i++)
					data[j][k][gridResolution-1-i]=newvol[k][j][i]; // true up-side-down

		System.out.println("Dimensions of the data structure: " + data.length + " " + data[0].length + " " + data[0][0].length);
    }

	// Create an octree structure to accelerate the rendering
	public void makeOctree()
	{
		int dimX = data[0][0].length;
		int dimY = data[0].length;
		int dimZ = data.length;
		int newDimX = (dimX / octreeBlockSize) + 1;
		int newDimY = (dimY / octreeBlockSize) + 1;
		int newDimZ = (dimZ / octreeBlockSize) + 1;
		octree = new int[dimY][dimX][newDimZ];

		for (int j = 0; j < newDimY; j++) 
			for (int i = 0; i < newDimX; i++)		
				for (int k = 0; k < newDimZ; k++) 
				{
					octree[j][i][k]=0;
					for (int j2 = 0; j2 < newDimY+1; j2++)
						for (int i2 = 0; i2 < newDimX+1; i2++)		
							for (int k2 = 0; k2 < newDimZ+1; k2++) 
							{
								int x = i * octreeBlockSize + i2;
								int y = j * octreeBlockSize + j2;
								int z = k * octreeBlockSize + k2;
								if (x < dimX && y < dimY && z < dimZ && octree[j][i][k] < data[z][y][x]) // if max value in the block
									octree[j][i][k] = data[z][y][x];
							}
				}
	}

	/**
	* This function returns the 3D gradient for the volumetric dataset (data variable). 
	* Note that the gradient values at the sides of the volume is not be computable. 
	* Each cell element containing a 3D vector, the result is therefore a 4D array.
	*/
	int [][][][] Gradient()
	{
		int[][][][] gradient=null;
		int dimX = data[0][0].length;
		int dimY = data[0].length;
		int dimZ = data.length;
		gradient = new int[dimZ-2][dimY-2][dimX-2][3]; //-2 due gradient not being computable at borders 
		for (int k = 1; k < dimZ-1; k++) 
			for (int j = 1; j < dimY-1; j++) 
				for (int i = 1; i < dimX-1; i++)
				{
						gradient[k-1][j-1][i-1][0] = (data[k][j][i+1] - data[k][j][i-1]) / 2;
						gradient[k-1][j-1][i-1][1] = (data[k][j+1][i] - data[k][j-1][i]) / 2;
						gradient[k-1][j-1][i-1][2] = (data[k+1][j][i] - data[k-1][j][i]) / 2;
				}
		return gradient;
	}

	// Interpolation functions
	// -----------------------
	// Trilinear interpolation functions
	double trilinearInterpolate(double x, double y, double z) {

        // Check if the position is within the grid boundaries
        if (x < 0 || y < 0 || z < 0 || x >= data[0][0].length-1 || y >= data[0].length-1 || z >= data.length-1) 
            return 0.0; // Return 0 if the position outside grid

		// Calculate the indices of the surrounding grid points
        int x0 = (int) Math.floor(x);
        int y0 = (int) Math.floor(y);
        int z0 = (int) Math.floor(z);
        // Calculate the fractional distances within the grid cell
        x = x - x0;
        y = y - y0;
        z = z - z0;
        // Perform trilinear interpolation
        double c =	data[z0][y0][x0] * (1-x) * (1-y) * (1-z) +
					data[z0][y0][x0+1] * (x) * (1-y) * (1-z) +
					data[z0][y0+1][x0] * (1-x) * (y) * (1-z) +
					data[z0][y0+1][x0+1] * (x) * (y) * (1-z) +
					data[z0+1][y0][x0] * (1-x) * (1-y) * (z) +
					data[z0+1][y0][x0+1] * (x) * (1-y) * (z) +
					data[z0+1][y0+1][x0] * (1-x) * (y) * (z) +
					data[z0+1][y0+1][x0+1] * (x) * (y) * (z);

        return c;
    }

	double[] trilinearInterpolateGrad(double x, double y, double z, int[][][][] gradient) {
	
		// Grad not computable at the borders
		x = x - 1;
		y = y - 1;
		z = z - 1;

		if  (x<=0 || y<=0 || z<=0 || x>=data[0][0].length-3 || y>=data[0].length-3 || z>=data.length-3 )
			return null; // Return null if position outside grid
		
		int x0 = (int) Math.floor(x);
		int y0 = (int) Math.floor(y);
		int z0 = (int) Math.floor(z);
		x = x - x0;
		y = y - y0;
		z = z - z0;
		
		double [] gradientInterpolation = new double[3];
		
		for(int i= 0; i < 3; i++) {

			gradientInterpolation[i] = gradient[z0][y0][x0][i] * (1-x) * (1-y) * (1-z) +
									   gradient[z0][y0][x0+1][i] * (x) * (1-y) * (1-z) +
									   gradient[z0][y0+1][x0][i] * (1-x) * (y) * (1-z) +
									   gradient[z0][y0+1][x0+1][i] * (x) * (y) * (1-z) +
									   gradient[z0+1][y0][x0][i] * (1-x) * (1-y) * (z) +
									   gradient[z0+1][y0][x0+1][i] * (x) * (1-y) * (z) +
									   gradient[z0+1][y0+1][x0][i] * (1-x) * (y) * (z) +
									   gradient[z0+1][y0+1][x0+1][i] * (x) * (y) * (z);
		}
		return gradientInterpolation;
	}
	// Cubic interpolation functions
	double cubicInterpolate(double x, double y, double z) {
		// Check if the position is within the grid boundaries
		if (x < 1 || y < 1 || z < 1 || x >= data[0][0].length - 2 || y >= data[0].length - 2 || z >= data.length - 2)
			return 0.0; // Return 0 if position outside the grid

		// Calculate the indices of the surrounding grid points
		int x0 = (int) Math.floor(x);
		int y0 = (int) Math.floor(y);
		int z0 = (int) Math.floor(z);

		// Calculate the fractional distances within the grid cell
		double dx = x - x0;
		double dy = y - y0;
		double dz = z - z0;

		// Perform cubic interpolation
		double[] values = new double[4];
		for (int i = 0; i < 4; i++) {
			values[i] = cubicInterpolate1D(data[z0 - 1 + i][y0 - 1][x0 - 1], data[z0 - 1 + i][y0][x0 - 1],
					data[z0 - 1 + i][y0 + 1][x0 - 1], data[z0 - 1 + i][y0 + 2][x0 - 1], dx);
		}
		return cubicInterpolate1D(values[0], values[1], values[2], values[3], dy) * (1 - dz) +
				cubicInterpolate1D(values[0], values[1], values[2], values[3], dy + 1) * dz;
	}

	double cubicInterpolate1D(double p0, double p1, double p2, double p3, double t) {
		double t2 = t * t;
		double t3 = t2 * t;
		double a0 = p3 - p2 - p0 + p1;
		double a1 = p0 - p1 - a0;
		double a2 = p2 - p0;
		double a3 = p1;

		return a0 * t3 + a1 * t2 + a2 * t + a3;
	}

	double[] cubicInterpolateGrad(double x, double y, double z, int[][][][] gradient) {
		// Check if the position is within the grid boundaries
		if (x < 1 || y < 1 || z < 1 || x >= data[0][0].length - 2 || y >= data[0].length - 2 || z >= data.length - 2)
			return null; // Return null if the position is outside the grid

		// Calculate the indices of the surrounding grid points
		int x0 = (int) Math.floor(x);
		int y0 = (int) Math.floor(y);
		int z0 = (int) Math.floor(z);

		// Calculate the fractional distances within the grid cell
		double dx = x - x0;
		double dy = y - y0;
		double dz = z - z0;

		// Perform cubic interpolation for each component of the gradient
		double[] gradientInterpolation = new double[3];
		for (int i = 0; i < 3; i++) {
			double[] values = new double[4];
			for (int j = 0; j < 4; j++) {
				values[j] = cubicInterpolate1D(gradient[z0 - 1 + j][y0 - 1][x0 - 1][i], gradient[z0 - 1 + j][y0][x0 - 1][i],
						gradient[z0 - 1 + j][y0 + 1][x0 - 1][i], gradient[z0 - 1 + j][y0 + 2][x0 - 1][i], dx);
			}
			gradientInterpolation[i] = cubicInterpolate1D(values[0], values[1], values[2], values[3], dy) * (1 - dz) +
					cubicInterpolate1D(values[0], values[1], values[2], values[3], dy + 1) * dz;
		}
		return gradientInterpolation;
	}
	// -----------------------

	// Rendering functions
	// -----------------------
	public int[][] renderIso(int [][][][] gradient, int isovalue)
	{
		int image[][] = new int[imgResolution][imgResolution];
		for (int j = 0; j < imgResolution; j++)
			for (int i = 0; i < imgResolution; i++)
			{
				image[j][i] = 0;
				double k = 0;
				while (k < data.length-1) 
				{
					// Keep going along the ray until we reach the isosurface
					// Less memory efficient than working by slices
					if (
						octree[((int) (j/zoom))/octreeBlockSize][((int) (i/zoom))/octreeBlockSize][(int) k/octreeBlockSize]>=isovalue
					)
					{
						double inter = trilinearInterpolate(i/zoom, j/zoom, k);
						if (inter > isovalue)
						{
							double g[] = trilinearInterpolateGrad(i/zoom, j/zoom, k, gradient);
							if (g != null)
							{
								//Normalise gradient before shading
								double norm = Math.sqrt(g[0] * g[0] + 
														g[1] * g[1] +
														g[2] * g[2]);
								if (norm > 0)
									image[j][i] = Math.min((int) Math.abs(255. * g[2] / norm), 255);
							}
						}
						k += samplingDistance;
					}
					else
						k += octreeBlockSize;
				}
			}
		return image;
	}

	public int[][] renderIsoWithoutOctree(int[][][][] gradient, int isovalue) {
		int[][] image = new int[imgResolution][imgResolution];
		for (int j = 0; j < imgResolution; j++) {
			for (int i = 0; i < imgResolution; i++) {
				image[j][i] = 0;
				double k = 0;
				while (k < data.length - 1) {
					if (trilinearInterpolate(i/zoom, j/zoom, k) > isovalue) {
						// Multiply zoom by 2 in trilinearInterpolateGrad()
						// (I dont know why, but tiral & error showed that it does not produce horizontal artifacts)
						// Though the shading is not perfect in this case
						double g[] = trilinearInterpolateGrad(i/zoom, j/zoom, k, gradient); 
						if (g != null)
						{
							//Normalise gradient before shading
							double norm = Math.sqrt(g[0] * g[0] + 
													g[1] * g[1] +
													g[2] * g[2]);
							if (norm > 0)
								image[j][i] = Math.min((int) Math.abs(255. * g[2] / norm), 255);
						}
					}
					k += samplingDistance;
				}
			}
		}
		return image;
	}

	int[][] ExtractSlice(int z)
	{
		return data[z];
	}

	int[][][] ExtractSliceOfGradient(int [][][][] gradient, int z, double amplification)
	{
		int slice[][][]=new int [gradient[0].length][gradient[0][0].length][3];
		for (int j=0;j<gradient[0].length;j++)
			for (int i=0;i<gradient[0][0].length;i++)
			{
				slice[j][i][0]=Math.min((int) Math.abs(gradient[z][j][i][0]*amplification),255);
				slice[j][i][1]=Math.min((int) Math.abs(gradient[z][j][i][1]*amplification),255);
				slice[j][i][2]=Math.min((int) Math.abs(gradient[z][j][i][2]*amplification),255);
			}
		return slice;
	}

	public int[][] Render(int [][][][] gradient, int isovalue, boolean positiveDirection)
	{
		int dimX = data[0][0].length;
		int dimY = data[0].length;
		int dimZ = data.length;
		int[][] image = new int[dimY][dimX];
		//Initialising the image values to 0 (black)
		for (int j = 0; j < dimY; j++) 
			for (int i = 0; i < dimX; i++)
				image[j][i]=0;
		for (int k = 0; k < dimZ-2; k++) 
		{
			//The next lines of code define the direction of projection.
			int sliceId; 
			if (positiveDirection)
				sliceId=k;
			else sliceId=dimZ-3-k;

			int[][] slice=ExtractSlice(sliceId);
			int[][][] gradientSlice=ExtractSliceOfGradient(gradient, sliceId,10.0);

			for (int j = 0; j < dimY-2; j++) 
				for (int i = 0; i < dimX-2; i++){
					if (slice[j][i]>isovalue)
					{
						//Normalising gradient before shading
						double norm=Math.sqrt(
							gradientSlice[j][i][0]*gradientSlice[j][i][0]+
							gradientSlice[j][i][1]*gradientSlice[j][i][1]+
							gradientSlice[j][i][2]*gradientSlice[j][i][2]);
						image[j][i]=Math.min((int) Math.abs(255.*gradientSlice[j][i][2]/norm),255);
					}
				}
		}
		return image;
	}
	// -----------------------

	// Swap the z axis with the x or y axis
	void SwapZAxis(int axis)
	{
		if (axis==2)
			return;
		int dimX=data[0][0].length;
		int dimY=data[0].length;
		int dimZ=data.length;
		int newvol[][][];
		if (axis==0)
		{
			newvol=new int[dimX][dimY][dimZ];
			for (int k = 0; k < dimZ; k++) 
				for (int j = 0; j < dimY; j++) 
					for (int i = 0; i < dimX; i++)
						newvol[i][j][k]=data[k][j][i];
		}
		else
		{
			newvol=new int[dimY][dimZ][dimX];
			for (int k = 0; k < dimZ; k++) 
				for (int j = 0; j < dimY; j++) 
					for (int i = 0; i < dimX; i++)
						newvol[j][k][i]=data[k][j][i];
		}
		data=newvol;
	}

}


public class Heart
{
    // Saving the image to a file
    public static void SaveImage(String filename, int[][] im) {
		BufferedImage image = new BufferedImage(im.length, im[0].length, BufferedImage.TYPE_BYTE_GRAY );
		for (int j = 0; j < im.length; j++) 
			for (int i = 0; i < im[0].length; i++) 
				image.setRGB(j, i, im[j][i]*256*256+im[j][i]*256+im[j][i]);
		
		File f = new File(filename);
		try 
		{
			ImageIO.write(image, "tiff", f);
			System.out.println("Image saved as "+filename);
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		}
    }


    // Main function
    public static void main(String[] args) {

        Volume v = new Volume();

        if (args.length < 3) {
            System.out.println("\n\nUsage: java Heart <grid_resolution> <image_resolution> <sampling_distance>\n\n");
            return;
        }

		System.out.println("\n\nVisualize Isosurface of Heart Equation!");
		
		v.setGridResolution(Integer.parseInt(args[0]));
        v.setImageResolution(Integer.parseInt(args[1]));
		v.setSamplingDistance(Float.parseFloat(args[2]));
		v.setZoom(); // Will be used to scale the volume to the image resolution

		System.out.println("\nSampling equation...");
		v.makeHeart(); // Sample the heart equation
		
		v.SwapZAxis(2); // No need, I've already done it in the makeHeart function
		
		v.makeOctree(); // For faster rendering, precompute the octree
		
		System.out.println("Computing gradient...");
        int[][][][] gradient = v.Gradient(); 
		
		System.out.println("Rendering isosurface...");
		System.out.println("Image resolution: " + v.getImageResolution() + "x" + v.getImageResolution() + " pixels");
		System.out.println("Sampling distance: " + v.getSamplingDistance()); 
		int [][] image = v.renderIso(gradient, 0); // Ray casting with octree
		//int [][] image = v.Render(gradient, 0, false); // Ray casting using 2D slices of the volume
		//int [][] image = v.renderIsoWithoutOctree(gradient, 0); // Ray casting without octree

		System.out.println("Saving image...");
		SaveImage("heart_isotest.tiff", image);
    }
}