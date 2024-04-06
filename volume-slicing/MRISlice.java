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
			System.out.println(fileName);
			FileInputStream f = new FileInputStream(fileName);
			DataInputStream d = new DataInputStream(f);

			d.readFully(dataBytes);
			
			//Copying the byte values into the floating-point array

			for (int k=0;k<sizeZ;k++)
				for (int j=0;j<sizeY;j++)
					for (int i=0;i<sizeX;i++)
						data[k][j][i]=dataBytes[k*sizeX*sizeY+j*sizeX+i+headerSize] & 0xff;
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
}

public class MRISlice
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
			System.out.println("Saving image");
			ImageIO.write(image, "png", f);
		}
		catch (IOException e) 
		{
			e.printStackTrace();
		}
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
			System.out.println("Saving image RGB\n");
			ImageIO.write(image, "png", f);
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		}
	}

	public static void main(String[] args) 
	{
		Volume v = new Volume();
		v.ReadData("./bighead_den256X256X225B62H.raw", Integer.parseInt(args[0]), Integer.parseInt(args[1]), Integer.parseInt(args[2]), Integer.parseInt(args[3]));
		int gradient[][][][] = v.Gradient();
		if (gradient==null)
			System.out.println("err\n");
		SaveImage("im.tiff",v.ExtractSlice(Integer.parseInt(args[4])));
		SaveImage("imGradientMagnitude.tiff",v.ExtractSliceOfGradientMagnitude(gradient,Integer.parseInt(args[4])));
		SaveImageRGB("imGradient.tiff",v.ExtractSliceOfGradient(gradient,Integer.parseInt(args[4])));
	}
}