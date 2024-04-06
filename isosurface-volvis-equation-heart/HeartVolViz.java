// import java.io.DataInputStream;
// import java.io.FileInputStream;
// import java.awt.Graphics;
import java.awt.*;
// import java.awt.Image;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;
import java.lang.Math;
// import java.awt.Color;

class Volume
{
    // Data structure to store the heart equation values
	double heartData[][][];

    // Helper function to evaluate the heart equation at a given point
    double evaluateEquation(double x, double y, double z) {
        return - (Math.pow((Math.pow(x, 2) + 2 * Math.pow(y, 2) + Math.pow(z, 2) - 1), 3) -
               Math.pow(x, 2) * Math.pow(z, 3) - 0.1 * Math.pow(y, 2) * Math.pow(z, 3));
    }

    // Helper function to map a value from one range to another
    double map(double value, double start1, double stop1, double start2, double stop2) {
        return start2 + (stop2 - start2) * ((value - start1) / (stop1 - start1));
    }

    // Task 1: Sampling the heart equation and storing it in a 3D array
    void sampleHeartEquation(int gridResolution) {

        heartData = new double[gridResolution][gridResolution][gridResolution];
        // Sample the heart equation within the grid
        for (int x = 0; x < gridResolution; x++) {
            double xPos = map(x, 0, gridResolution - 1, -2, 2);
            for (int y = 0; y < gridResolution; y++) {
                double yPos = map(y, 0, gridResolution - 1, -2, 2);
                for (int z = 0; z < gridResolution; z++) {
                    double zPos = map(z, 0, gridResolution - 1, -2, 2);
                    double equationResult = evaluateEquation(xPos, yPos, zPos);
                    heartData[x][y][z] = equationResult;
                }
            }
        }
    }

    // Task 2: Basic visualization using image-order algorithm
    BufferedImage visualizeBasic(int resolution) {
        System.out.println("Resolution: " + resolution);
        BufferedImage image = new BufferedImage(resolution, resolution, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g2d = image.createGraphics();

        // Render the isosurface
        for (int x = 0; x < resolution - 1; x++) {
            for (int y = 0; y < resolution - 1; y++) {
                for (int z = 0; z < resolution - 1; z++) {
                    double value1 = heartData[x][y][z];
                    double value2 = heartData[x + 1][y][z];
                    double value3 = heartData[x][y + 1][z];
                    double value4 = heartData[x + 1][y + 1][z];
                    double value5 = heartData[x][y][z + 1];
                    double value6 = heartData[x + 1][y][z + 1];
                    double value7 = heartData[x][y + 1][z + 1];
                    double value8 = heartData[x + 1][y + 1][z + 1];

                    // Check if the isosurface intersects the current cube
                    if ((value1 >= 0 && value2 < 0) || (value1 < 0 && value2 >= 0) ||
                            (value3 >= 0 && value4 < 0) || (value3 < 0 && value4 >= 0) ||
                            (value5 >= 0 && value6 < 0) || (value5 < 0 && value6 >= 0) ||
                            (value7 >= 0 && value8 < 0) || (value7 < 0 && value8 >= 0)) {
                        // Interpolate the position of the isosurface within the cube
                        // double t = interpolate(value1, value2);
                        double t = trilinearInterpolation(x, y, z, resolution);
                        double xInterpolated = map(x + t, 0, resolution - 1, -2, 2);
                        // double yInterpolated = map(y, 0, resolution - 1, -2, 2);
                        double zInterpolated = map(z, 0, resolution - 1, -2, 2);

                        // Map the interpolated position to the image coordinates
                        int pixelX = (int) map(xInterpolated, -2, 2, 0, resolution - 1);
                        // int pixelY = (int) map(yInterpolated, -2, 2, 0, resolution - 1);
                        int pixelZ = (int) map(zInterpolated, -2, 2, 0, resolution - 1);

                        // Set the pixel color to visualize the isosurface
                        g2d.setColor(Color.GRAY);
                        g2d.fillRect(pixelX, pixelZ, 1, 1);
                    }
                }
            }
        }
        return image;
    }

    // Function to create a basic visualization of the volume shape
    BufferedImage visualizeVolume(int resolution) {
        BufferedImage image = new BufferedImage(resolution, resolution, BufferedImage.TYPE_INT_RGB);

        // Find the maximum value in the volume data
        double maxVal = Double.MIN_VALUE;
        for (int i = 0; i < resolution; i++) {
            for (int j = 0; j < resolution; j++) {
                for (int k = 0; k < resolution; k++) {
                    maxVal = Math.max(maxVal, heartData[i][j][k]);
                }
            }
        }

        // Create the visualization
        for (int i = 0; i < resolution; i++) {
            for (int j = 0; j < resolution; j++) {
                for (int k = 0; k < resolution; k++) {
                    int intensity = (int) (255 * (heartData[i][j][k] / maxVal)); // Scale the intensity
                    int rgb = (intensity << 16) | (intensity << 8) | intensity; // Grayscale color
                    image.setRGB(i, j, rgb);
                }
            }
        }
        return image;
    }


    // Task 3: Shading using gradient vectors
    double[] calculateGradient(int x, int y, int z, int resolution) {
        double dx = (heartData[Math.min(x + 1, resolution - 1)][y][z] - heartData[Math.max(x - 1, 0)][y][z]) / 2.0;
        double dy = (heartData[x][Math.min(y + 1, resolution - 1)][z] - heartData[x][Math.max(y - 1, 0)][z]) / 2.0;
        double dz = (heartData[x][y][Math.min(z + 1, resolution - 1)] - heartData[x][y][Math.max(z - 1, 0)]) / 2.0;
        double magnitude = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (magnitude != 0) {
            dx /= magnitude;
            dy /= magnitude;
            dz /= magnitude;
        }
        return new double[]{dx, dy, dz};
    }

    // Task 4: Trilinear interpolation for image quality improvement
    double trilinearInterpolation(double x, double y, double z, int resolution) {

        // Calculate the indices of the surrounding grid points
        int x0 = (int) Math.floor(x);
        int x1 = x0 + 1;
        int y0 = (int) Math.floor(y);
        int y1 = y0 + 1;
        int z0 = (int) Math.floor(z);
        int z1 = z0 + 1;

        // Check if the position is within the grid boundaries
        if (x0 < 0 || x1 >= resolution || y0 < 0 || y1 >= resolution || z0 < 0 || z1 >= resolution) {
            return 0.0; // Return 0 if the position is outside the grid
        }

        // Calculate the fractional distances within the grid cell
        double xd = x - x0;
        double yd = y - y0;
        double zd = z - z0;

        // Perform trilinear interpolation
        double c00 = heartData[x0][y0][z0] * (1 - xd) + heartData[x1][y0][z0] * xd;
        double c01 = heartData[x0][y0][z1] * (1 - xd) + heartData[x1][y0][z1] * xd;
        double c10 = heartData[x0][y1][z0] * (1 - xd) + heartData[x1][y1][z0] * xd;
        double c11 = heartData[x0][y1][z1] * (1 - xd) + heartData[x1][y1][z1] * xd;
        double c0 = c00 * (1 - yd) + c10 * yd;
        double c1 = c01 * (1 - yd) + c11 * yd;
        double c = c0 * (1 - zd) + c1 * zd;

        return c;
    }

    // Task: Image visualization using shading
    BufferedImage visualizeShading(int resolution) {
        BufferedImage image = new BufferedImage(resolution, resolution, BufferedImage.TYPE_INT_ARGB);
        for (int i = 0; i < resolution; i++) {
            double x = -2 + i * 4.0 / resolution;
            for (int j = 0; j < resolution; j++) {
                double y = -2 + j * 4.0 / resolution;
                for (int k = 0; k < resolution; k++) {
                    double z = -2 + k * 4.0 / resolution;
                    double value = trilinearInterpolation(x, y, z, resolution);
                    if (value >= 0) {
                        double[] gradient = calculateGradient(i, j, k, resolution);
                        double shading = Math.abs(gradient[0]) + Math.abs(gradient[1]) + Math.abs(gradient[2]);
                        if (shading < 2) {
                            int gray = (int) (255 * (1 - value));
                            int rgb = (gray << 16) | (gray << 8) | gray;
                            image.setRGB(i, j, rgb);
                        }
                    }
                }
            }
        }
        return image;
    }
	
	/**
	* This function returns the 3D gradient for the volumetric dataset (data variable). Note that the gradient values at the sides of the volume is not be computable. Each cell element containing a 3D vector, the result is therefore a 4D array.
	*/
	double [][][][] Gradient()
	{
		double[][][][] gradient=null;
		int dimX=heartData[0][0].length;
		int dimY=heartData[0].length;
		int dimZ=heartData.length;
		gradient=new double[dimZ-2][dimY-2][dimX-2][3]; //-2 due gradient not being computable at borders 
		for (int k = 1; k < dimZ-1; k++) 
			for (int j = 1; j < dimY-1; j++) 
				for (int i = 1; i < dimX-1; i++)
				{
						gradient[k-1][j-1][i-1][0]=(heartData[k][j][i+1]-heartData[k][j][i-1])/2;
						gradient[k-1][j-1][i-1][1]=(heartData[k][j+1][i]-heartData[k][j-1][i])/2;
						gradient[k-1][j-1][i-1][2]=(heartData[k+1][j][i]-heartData[k-1][j][i])/2;
				}
		return gradient;
	}

	public int[][] Render(int [][][][] gradient) 
    {
		int dimX=gradient[0][0].length;
        int dimY=gradient[0].length;
        int dimZ=gradient.length;
        int[][] image=new int[dimZ][dimY];
        for (int k = 0; k < dimZ; k++) 
            for (int j = 0; j < dimY; j++) 
                for (int i = 0; i < dimX; i++)
                {
                    double shading = Math.abs(gradient[k][j][i][0]) + Math.abs(gradient[k][j][i][1]) + Math.abs(gradient[k][j][i][2]);
                    if (shading < 2) {
                        int gray = (int) (255 * (1 - heartData[k][j][i]));
                        image[k][j]=gray;
                    }
                }
        return image;
	}

}

public class HeartVolViz
{
    // Saving the image to a file
    public static void SaveImage(String filename, BufferedImage image) {
        File file = new File(filename);
        try {
            ImageIO.write(image, "png", file);
            System.out.println("Image saved successfully as " + filename);
        } catch (IOException e) {
            System.out.println("An error occurred while saving the image: " + e.getMessage());
            e.printStackTrace();
        }
    }


    // Main function
    public static void main(String[] args) {
        Volume v = new Volume();
        if (args.length < 2) {
            System.out.println("Usage: java CW <grid_resolution> <image_resolution>");
            return;
        }
        int gridResolution = Integer.parseInt(args[0]);
        int imageResolution = Integer.parseInt(args[1]);

        v.sampleHeartEquation(gridResolution);
        BufferedImage image = v.visualizeBasic(imageResolution);
        BufferedImage image2 = v.visualizeVolume(imageResolution);
        BufferedImage image3 = v.visualizeShading(imageResolution);
        SaveImage("result.tiff", image);
        SaveImage("result2.tiff", image2);
        SaveImage("result3.tiff", image3);
    }
}