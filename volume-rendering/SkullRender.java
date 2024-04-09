import java.io.DataInputStream;
import java.io.FileInputStream;
// import java.awt.Graphics;
// import java.awt.Image;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;
import java.lang.Math;
import java.awt.Color;

class Volume
{
	int data[][][];
	float zoom = 2F;
	int resolution = 512;
	
	public int GetResolution() 
	{
		return resolution;
	}
	
	public void SetResolution(int res) 
	{
		resolution = res;
	}
	
	public float GetZoom() 
	{
		return zoom;
	}
	
	public void SetZoom(float z) 
	{
		zoom = z;
	}
	
	/**
	* This function reads a volume dataset from disk and put the result in the data array
	* @param amplification allows increasing the brightness of the slice by a constant.
	*/
	boolean ReadData(String fileName, int sizeX, int sizeY, int sizeZ, int headerSize)
	{
		int cpt=0;
		byte dataBytes[]=new byte [sizeX*sizeY*sizeZ+headerSize];
		data = new int[sizeZ][sizeY][sizeX];
	    try
		{
			FileInputStream f = new FileInputStream(fileName);
			DataInputStream d = new DataInputStream(f);

			d.readFully(dataBytes);
			
			//Copying the byte values into the floating-point array

			for (int k=0;k<sizeZ;k++)
				for (int j=0;j<sizeY;j++)
					for (int i=0;i<sizeX;i++)
						data[k][j][i]=dataBytes[k*256*sizeY+j*sizeX+i+headerSize] & 0xff;
			d.close();
		}
		catch(Exception e)
		{ 
			System.out.println("Exception : "+cpt+e);
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
	* This function returns the 3D gradient vector for the zth slice of the volume gradient. Values may be positive-only and capped to 255. 
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
	* This function returns the 3D gradient for the volumetric dataset (data variable). Note that the gradient values at the sides of the volume is not be computable. Each cell element containing a 3D vector, the result is therefore a 4D array.
	*/
	int [][][][] Gradient()
	{
		int[][][][] gradient=null;
		int dimX=data[0][0].length;
		int dimY=data[0].length;
		int dimZ=data.length;
		gradient=new int[dimZ-2][dimY-2][dimX-2][3]; //-2 due gradient not being computable at borders 
		for (int k = 1; k < dimZ-1; k++) 
			for (int j = 1; j < dimY-1; j++) 
				for (int i = 1; i < dimX-1; i++)
				{
						gradient[k-1][j-1][i-1][0]=(data[k][j][i+1]-data[k][j][i-1])/2;
						gradient[k-1][j-1][i-1][1]=(data[k][j+1][i]-data[k][j-1][i])/2;
						gradient[k-1][j-1][i-1][2]=(data[k+1][j][i]-data[k-1][j][i])/2;
				}
		return gradient;
	}

	double computeTrilinearIntepolation(double x0, double y0,double z0) {
		
		int x = (int) x0;
        int y = (int) y0;
        int z = (int) z0;
		double xd = x0 - x;
        double yd = y0 - y;
        double zd = z0 - z;
		
		double c1 = data[z][y][x]*(1-zd)*(1-yd)*(1-xd)+
					data[z+1][y][x]*zd*(1-yd)*(1-xd)+ 
					data[z][y+1][x]*(1-zd)*yd*(1-xd)+
					data[z+1][y+1][x]*zd*yd*(1-xd);
		
		double c2 = data[z][y][x+1]*(1-zd)*(1-yd)*xd+
					data[z+1][y][x+1]*zd*(1-yd)*xd+
					data[z][y+1][x+1]*(1-zd)*yd*xd+
					data[z+1][y+1][x+1]*zd*yd*xd;
		
		return c1+c2;
	}

	double[] computeTrilinearIntepolationGradient(double x, double y, double z, int[][][] gradientSlice) {
	
		// Get the integer coordinates of the 8 surrounding voxels
		int x0 = (int)x;
		int x1 = x0 + 1;
		int y0 = (int)y;
		int y1 = y0 + 1;
		int z0 = (int)z;
		int z1 = z0 + 1;
	
	
		double xd = x - x0;
		double yd = y - y0;
		double zd = z - z0;
		
		double [] gradientInterpolation = new double[3];
		
		for(int i= 0; i<gradientInterpolation.length;i++) {
			double c00 = gradientSlice[y0][x0][i] + gradientSlice[y0][x0][i] * xd;
			
			double c01 = gradientSlice[y1][x0][i] * (1 - xd) + gradientSlice[y1][x0][i] * xd;
			double c10 = gradientSlice[y0][x1][i] * (1 - xd) + gradientSlice[y0][x1][i] * xd;
			double c11 = gradientSlice[y1][x1][i] * (1 - xd) + gradientSlice[y1][x1][i] * xd;
			double c0 = c00 * (1 - yd) + c10 * yd;
			double c1 = c01 * (1 - yd) + c11 * yd;
			double c = c0 * (1 - zd) + c1 * zd;
	
			gradientInterpolation[i] =   c;
		}
	
		
	
		return gradientInterpolation;
	}


	/**
	* This function returns an image of a isosurface visualisation projected along the z axis.
	* @param direction The direction of the ray along the axis
	* @param isovalue The threshold value for delimitating the isosurface
	*/
	public int[][] Render(int [][][][] gradient, int isovalue, boolean positiveDirection)
	{
		//The algorithm will work by projecting slices along the z-axis, re-using code from the previous lab
		//This also has  the property to work in memory order.
		int res = GetResolution();
		float zoom = GetZoom();
		int dimX=data[0][0].length-1;
		int dimY=data[0].length-1;
		int dimZ=data.length-1;
		// System.out.println(dimX);
		// System.out.println(dimY);
		// System.out.println(dimZ);
		
		int[][] image = new int[res][res];
		for (int j = 0; j < res-2; j++) 
			for (int i = 0; i < res-2; i++)
				image[j][i]=0;		
		double z= 0;
		while (z< dimZ-3)
		
		{
			//The next lines of code define the direction of projection.
			int sliceId; 
			int k= (int) z;
			if (positiveDirection)
				sliceId=k;
			else sliceId=dimZ-3-k;
					
			int[][][] gradientSlice=ExtractSliceOfGradient(gradient,sliceId,10.0);

			double y= 0;
					
			while (y< dimY-3)
						
			{
				double x=0;
				while(x<dimX-3)
							
				{
					int i= (int) (x * zoom);
					int j= (int) (y * zoom);
					double sample= computeTrilinearIntepolation(x,y,z);
					if (sample>isovalue)
					{	
						double [] interpolatedGradient = computeTrilinearIntepolationGradient(x,y,z,gradientSlice);
						double norm = Math.sqrt(interpolatedGradient[0]*interpolatedGradient[0]+interpolatedGradient[1]*interpolatedGradient[1]+interpolatedGradient[2]*interpolatedGradient[2]);
						//Normalise gradient before shading
						// double norm=Math.sqrt(gradientSlice[j][i][0]*gradientSlice[j][i][0]+gradientSlice[j][i][1]*gradientSlice[j][i][1]+gradientSlice[j][i][2]*gradientSlice[j][i][2]);
						image[j][i]=Math.min((int) Math.abs(255.*interpolatedGradient[2]/norm),255);
						//image[j][i]=Math.min((int) Math.abs(255.*gradientSlice[j][i][0]),255);
					}
					x+=(1/zoom);	
				}
				y+=(1/zoom);	
			}
			z+=0.5;
		}
		return image;
	}

	/**
	* This function swaps the x or y dimension with the z one, allowing projection on other faces of the volume.
	*/
	void SwapDimensions(int axis)
	{
		// TO-DO For projection from left of the volume
		// This function is broken
		if (axis == -1)
		{
			// Reverse array elements in the 3D array
			int layers = data.length;
			int rows = data[0].length;
			int columns = data[0][0].length;
			System.out.println("Layers: " + layers + " Rows: " + rows + " Columns: " + columns + "\n");
			for (int k = 0; k < layers; k++) {
				for (int j = 0; j < rows / 2; j++) {
					for (int i = 0; i < columns; i++) {
						// Swap elements along each dimension
						int temp = data[j][i][k];
						System.out.println(rows - j);
						data[j][i][k] = data[rows - j - 1][i][k];
						data[rows - j - 1][i][k] = temp;
					}
				}
			}  
			axis = 1;
		}
		if (axis == 0)
		{
			int[][][] newData = new int[data[0].length][data[0][0].length][data.length];
			for (int k = 0; k < data.length; k++) 
				for (int i = 0; i < data[0].length; i++) 
					for (int j = 0; j < data[0][0].length; j++)
						newData[j][i][k]=data[k][i][j];
			data=newData;
		}
		else if (axis == 1)
		{
			int[][][] newData = new int[data[0].length][data[0][0].length][data.length];
			for (int k = 0; k < data.length; k++) 
				for (int j = 0; j < data[0].length; j++) 
					for (int i = 0; i < data[0][0].length; i++)
						newData[j][i][k]=data[j][k][i];
			data=newData;
		}
		else
			return;
	}

}

public class SkullRender
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
		// A command line example: java SkullRender 256 256 225 62 95 0 false
		Volume v=new Volume();
		v.ReadData("./bighead_den256X256X225B62H.raw",Integer.parseInt(args[0]),Integer.parseInt(args[1]),Integer.parseInt(args[2]),Integer.parseInt(args[3]));
		v.SwapDimensions(Integer.parseInt(args[5]));
		int[][][][] gradient = v.Gradient();
		SaveImage("skull_iso.tiff",v.Render(gradient,Integer.parseInt(args[4]),Boolean.parseBoolean(args[6])));
	}
}