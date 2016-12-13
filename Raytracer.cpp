#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

void writePPM(const char * imagepath, int width, int height, int*** img, int color = 255);

int main()
{
	int width = 256*2, height = 256*2;
	double freq = 3;
	double thickness = 0.2;
	double f = 20;
	double proj_x = 3;
	double e_x = 25, e_y = 50, e_z = 25;
	int p1_r = 128, p1_g = 0, p1_b = 0;
	int p2_r = 0, p2_g = 128, p2_b = 0;
	int p3_r = 0, p3_g = 0, p3_b = 128;
	int kd = 40, ka = 10;
	// plane normal -- normalized
	double angle = 45;
	double n_x = cos(M_PI * angle / 180), n_y = 0, n_z = sin(M_PI * angle / 180), D = -f * sin(M_PI * angle / 180);
	//double pn3_x = cos(M_PI * angle / 180), pn3_y = 0, pn3_z = -sin(M_PI * angle / 180), p3_D = 35 * sin(M_PI * angle / 180);
	//double pn2_x = -cos(M_PI * angle / 180), pn2_y = 0, pn2_z = -sin(M_PI * angle / 180), p2_D = 35 * sin(M_PI * angle / 180);
	double pn3_x = 0, pn3_y = 0, pn3_z = 1, p3_D = 0;
	double pn2_x = 1, pn2_y = 0, pn2_z = 0, p2_D = 0;
	double pn1_x = 0, pn1_y = 1, pn1_z = 0, p1_D = 0;

	int*** arr = new int**[width];
	for (int i = 0; i < width; i++)
	{
		arr[i] = new int*[height];
		for (int j = 0; j < height; j++)
			arr[i][j] = new int[3]();
	}
	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++)
		{
			double x = (((j - 0)*(f - proj_x - proj_x)) / (height - 1)) + proj_x;
			double y = height*0.75-i;
			double dd_x = x - e_x, dd_y = y - e_y;
			double z = -1*(x*n_x + y*n_y + D)/n_z;
			double dd_z = z - e_z;
			double norm = sqrt(dd_x*dd_x + dd_y*dd_y + dd_z*dd_z);
			double d_x = dd_x / norm;
			double d_y = dd_y / norm;
			double d_z = dd_z / norm;
			double t1;
			double denom1 = (pn1_x*d_x + pn1_y*d_y + pn1_z*d_z);
			if (abs(denom1) - 1e-5 <= 0)
				t1 = 1e7;
			t1 = -(p1_D + (pn1_x*e_x + pn1_y*e_y + pn1_z*e_z)) / denom1;
			double t2;
			double denom2 = (pn2_x*d_x + pn2_y*d_y + pn2_z*d_z);
			if (abs(denom2) - 1e-5 <= 0)
				t2 = 1e7;
			t2 = -(p2_D + (pn2_x*e_x + pn2_y*e_y + pn2_z*e_z)) / denom2;
			double t3;
			double denom3 = (pn3_x*d_x + pn3_y*d_y + pn3_z*d_z);
			if (abs(denom3) - 1e-5 <= 0)
				t3 = 1e7;
			t3 = -(p3_D + (pn3_x*e_x + pn3_y*e_y + pn3_z*e_z)) / denom3;
			double denom = 0;
			double t = 0;
			int p_r, p_g, p_b;
			if (t1 >= 0 && (t1 <= t2 || t2 <= 0) && (t1 <= t3 || t3 <= 0))
			{
				denom = denom1;
				double intersect_x = e_x + t1*d_x;
				double intersect_y = e_y + t1*d_y;
				double intersect_z = e_z + t1*d_z;
				if (abs(cos(freq / f*intersect_x * 2 * M_PI)) < thickness || abs(cos(freq / f*intersect_z * 2 * M_PI)) < thickness)
				{
					p_r = 0;
					p_g = 0;
					p_b = 0;
				}
				else
				{
					p_r = 255;
					p_g = 255;
					p_b = 255;
				}
				t = t1;
			}
			else if (t2 >= 0 && (t2 <= t1 || t1 <= 0) && (t2 <= t3 || t3 <= 0))
			{
				//cout << "t2";
				denom = denom2;
				double intersect_x = e_x + t2*d_x;
				double intersect_y = e_y + t2*d_y;
				double intersect_z = e_z + t2*d_z;
				if (abs(cos(freq / f*intersect_z * 2 * M_PI)) < thickness || abs(cos(freq / height*intersect_y * 2 * M_PI)) < thickness)
				{
					p_r = 0;
					p_g = 0;
					p_b = 0;
				}
				else
				{
					p_r = 255;
					p_g = 255;
					p_b = 255;
				}
				t = t2;
			}
			else if (t3 >= 0 && (t3 <= t1 || t1 <= 0) && (t3 <= t2 || t2 <= 0))
			{
				//cout << "t3";
				denom = denom3;
				double intersect_x = e_x + t3*d_x;
				double intersect_y = e_y + t3*d_y;
				double intersect_z = e_z + t3*d_z;
				if (abs(cos(freq / f*intersect_x * 2 * M_PI)) < thickness || abs(cos(freq / height * intersect_y * 2 * M_PI)) < thickness)
				{
					p_r = 0;
					p_g = 0;
					p_b = 0;
				}
				else
				{
					p_r = 255;
					p_g = 255;
					p_b = 255;
				}
				/*
				p_r = p3_r;
				p_g = p3_g;
				p_b = p3_b;
				*/
				t = t3;
			}
			else
				continue;
			arr[i][j][0] = ka + kd*abs(denom)*p_r;
			arr[i][j][1] = ka + kd*abs(denom)*p_g;
			arr[i][j][2] = ka + kd*abs(denom)*p_b;

			// clip values
			arr[i][j][0] = arr[i][j][0] > 255 ? 255 : arr[i][j][0];
			arr[i][j][1] = arr[i][j][1] > 255 ? 255 : arr[i][j][1];
			arr[i][j][2] = arr[i][j][2] > 255 ? 255 : arr[i][j][2];
		}
	writePPM("test.ppm", width, height, arr);
	return 0;
}

void writePPM(const char * imagepath, int width, int height, int*** img, int color)
{
	std::ofstream f(imagepath);
	// header
	f << "P3" << "\n";
	// meta
	f << width << " " << height << " " << color << "\n";
	// data
	int r = 255, g = 0, b = 0;
	for (int i = 0; i < width; i++)
		for (int j = 0; j < height; j++)
			f << img[i][j][0] << " " << img[i][j][1] << " " << img[i][j][2] << "\n";

	f.close();
}
