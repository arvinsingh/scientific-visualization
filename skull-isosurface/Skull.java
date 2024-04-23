import java.io.DataInputStream;
import java.io.FileInputStream;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;
import java.lang.Math;
import java.awt.Color;
// import java.awt.Graphics;
// import java.awt.Image;

class Volume
{
	int data[][][];
	int imgResolution = 1024;
    float zoom = 4.f;
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

	public float getZoom() {
		return zoom;
	}

	public void setZoom(int z) {
		this.zoom = z;
	}

	public float getSamplingDistance() {
		return samplingDistance;
	}

	public void setSamplingDistance(float s) {
		this.samplingDistance = s;
	}

	/**
	* This function reads a volume dataset from disk and put the result in the data array
	* @param amplification allows increasing the brightness of the slice by a constant.
	*/
	boolean ReadData(String fileName, int sizeX, int sizeY, int sizeZ, int headerSize)
	{
		int cpt = 0;
		byte dataBytes[] = new byte [sizeX * sizeY * sizeZ + headerSize];
		data = new int[sizeZ][sizeY][sizeX];
	    try
		{
			FileInputStream f = new FileInputStream(fileName);
			DataInputStream d = new DataInputStream(f);

			d.readFully(dataBytes);
			
			//Copying the byte values into the floating-point array
			for (int k=0; k<sizeZ; k++)
				for (int j=0; j<sizeY; j++)
					for (int i=0; i<sizeX; i++)
						data[k][j][i] = dataBytes[k * 256 * sizeY + j * sizeX + i + headerSize] & 0xff;

			d.close();
		}
		catch(Exception e)
		{ 
			System.out.println("Exception : " + cpt + e);
			return false;
		}
		return true;
	}
	
	int[][] ExtractSlice(int z)
	{
		return data[z];
	}
	
	int[][] ExtractSliceOfGradientMagnitude(int [][][][] gradient, int z)
	{
		return ExtractSliceOfGradientMagnitude(gradient, z, 5.);
	}

	/**
	* This function returns the gradient magnitude for the zth slice of the volume gradient. 
	* @param amplification allows increasing the brightness of the slice by a constant.
	*/
	int[][] ExtractSliceOfGradientMagnitude(int [][][][] gradient, int z, double amplification)
	{
		int slice[][]=new int [gradient[0].length][gradient[0][0].length];
		for (int j=0;j<gradient[0].length;j++)
			for (int i=0;i<gradient[0][0].length;i++)
				slice[j][i]=(int) Math.min(Math.sqrt((gradient[z][j][i][0]*gradient[z][j][i][0]+gradient[z][j][i][1]*gradient[z][j][i][1]+gradient[z][j][i][2]*gradient[z][j][i][2]))*amplification,255); //min operator to make sure we are not having out of range pixel values
		return slice;
	}
	
	int[][][] ExtractSliceOfGradient(int [][][][] gradient, int z)
	{
		return 	ExtractSliceOfGradient(gradient, z, 15.0);
	}
	
	/**
	* This function returns the 3D gradient vector for the zth slice of the volume gradient. 
	* Values may be positive-only and capped to 255. 
	* @param amplification allows increasing the brightness of the slice by a constant.
	*/
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

	/**
	* This function returns an image of a isosurface visualisation projected along the z axis.
	* @param direction The direction of the ray along the axis
	* @param isovalue The threshold value for delimitating the isosurface
	*/
	public int[][] render(int [][][][] gradient, int isovalue, boolean positiveDirection)
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

	public int[][] renderContour(int [][][][] gradient, int isovalue)
	{
		int image[][] = new int[imgResolution][imgResolution];
		for (int j = 0; j < imgResolution; j++)
			for (int i = 0; i < imgResolution; i++)
			{
				double sum = 0;
				double k = 0;
				while (k < data.length-1) 
				{
					double inter = trilinearInterpolate(i/zoom, j/zoom, k);
					if (inter > isovalue)
					{	
						double g[] = trilinearInterpolateGrad(i/zoom, j/zoom, k, gradient);
						if (g != null)
						{
							double norm = g[0] * g[0] +
										  g[1] * g[1] + 
										  g[2] * g[2];
							sum += norm;
						}
					}
					k += samplingDistance;
				}
				image[j][i] = (int) (255.*sum / (20000. + sum)); //Scaling to range [0,255]
			}
		return image;
	}

	// Swap the z axis with the x or y axis
	void SwapZAxis(int axis)
	{
		if (axis == 2)
			return;
		int dimX = data[0][0].length;
		int dimY = data[0].length;
		int dimZ = data.length;
		int newvol[][][];
		if (axis == 0)
		{
			newvol = new int[dimX][dimY][dimZ];
			for (int k = 0; k < dimZ; k++) 
				for (int j = 0; j < dimY; j++) 
					for (int i = 0; i < dimX; i++)
						newvol[i][j][k]=data[k][j][i];
		}
		else
		{
			newvol = new int[dimY][dimZ][dimX];
			for (int k = 0; k < dimZ; k++) 
				for (int j = 0; j < dimY; j++) 
					for (int i = 0; i < dimX; i++)
						newvol[j][k][i]=data[k][j][i];
		}
		data = newvol;
	}

}

public class Skull
{
	public static void SaveImage(String name, int[][] im)
	{
		BufferedImage image = new BufferedImage(im.length, im[0].length, BufferedImage.TYPE_BYTE_GRAY );
		for (int j = 0; j < im.length; j++) 
			for (int i = 0; i < im[0].length; i++) 
				image.setRGB(j, i, im[j][i]*256*256+im[j][i]*256+im[j][i]);
		
		File f = new File(name);
		try 
		{
			ImageIO.write(image, "png", f);
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		}
		System.out.println("Image saved as "+name);
	}

	public static void SaveImageRGB(String name, int[][][] im)
	{
		BufferedImage image = new BufferedImage(im.length, im[0].length, BufferedImage.TYPE_INT_RGB );
		for (int j = 0; j < im.length; j++) 
			for (int i = 0; i < im[0].length; i++) 
			{
				Color c=new Color(Math.abs(im[j][i][0]),Math.abs(im[j][i][1]),Math.abs(im[j][i][2]));
				image.setRGB(j, i, c.getRGB());
			}
		
		File f = new File(name);
		try 
		{
			ImageIO.write(image, "png", f);
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		}
	}

	public static void main(String[] args) 
	{
		// Args: width, height, depth, header_size, isovalue, projection_axis, direction
		// A command line example: java Skull 95 0 false
		Volume v = new Volume();

		if (args.length < 3) {
            System.out.println("Usage: java Skull <isovalue> <projection_axis> <direction>");
            return;
        }

		System.out.println("Reading bighead_den256X256X225B62H.raw ...");
		System.out.println("Width x Height x Depth of the volume: 256x256x225");
		System.out.println("Header Size: 62 bytes");

		int imgWidth = 256;
		int imgHeight = 256;
		int imgDepth = 225;
		int headerSize = 62;

		v.ReadData("./bighead_den256X256X225B62H.raw", imgWidth, imgHeight, imgDepth, headerSize);
		v.SwapZAxis(Integer.parseInt(args[1]));
		v.makeOctree();

		System.out.println("Computing Gradient ...");
		int[][][][] gradient = v.Gradient();

		System.out.println("Isovalue: " + args[0]);
		System.out.println("Rendering Isosurface ...");
		int [][] image = v.renderIso(gradient, Integer.parseInt(args[0]));
		SaveImage("img/skull_iso_" + Integer.parseInt(args[0]) + ".png",image);
	
		System.out.println("Rendering Contour ...");
		int [][] contour = v.renderContour(gradient, Integer.parseInt(args[0]));
		SaveImage("img/skull_contour_" + Integer.parseInt(args[0]) + ".png", contour);

		System.out.println("Rendering using slices along the z axis ...");
		int [][] image2 = v.render(gradient, Integer.parseInt(args[0]), Boolean.parseBoolean(args[2]));
		SaveImage("skull_iso_using_slices_" + Integer.parseInt(args[0]) + ".png", image2);
	}
}