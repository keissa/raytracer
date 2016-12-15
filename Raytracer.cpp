#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>

#include "Vect.h"

#define OBJS_NUM 6

using namespace std;

void writePPM(const char * imagepath, int width, int height, int*** img, int color = 255);

bool intersect_plane(Vect pos, Vect ray_direction, Vect plane_normal, double D, two_doubles& td);
bool intersect_sphere(Vect pos, Vect ray_direction, Vect center, double radius, two_doubles& td);

int width = 256 * 2, height = 256 * 2;
double freq = 200;
double thickness = 0.2;
double f = 256;
double proj_x = 3;

Vect O(0, 0, 0);
Vect X(1, 0, 0);
Vect Y(0, 1, 0);
Vect Z(0, 0, 1);

Vect eye_pos = Vect(5+3, 3+0.5, 5+3);
Vect look_at = Vect(0, 0.5, 0);
// Camera three axes are cam_dir, cam_right, cam_down
Vect cam_dir = (look_at - eye_pos).normalize();
Vect cam_right = Y.crossProduct(cam_dir).normalize();
Vect cam_down = (cam_dir.crossProduct(cam_right).normalize())*(-1);

double max(double x, double y)
{
	return x >= y ? x : y;
}

two_doubles max(double x, double y, double z)
{
	two_doubles td;
	if (x >= y && x >= z)
	{
		td.denom = x;
		if (y >= z)
			td.t = y;
		else
			td.t = z;
	}
	if (y >= x && y >= z)
	{
		td.denom = y;
		if (x >= z)
			td.t = x;
		else
			td.t = z;
	}
	if (z >= y && z >= x)
	{
		td.denom = z;
		if (y >= x)
			td.t = y;
		else
			td.t = x;
	}
	return td;
}

int main()
{
	int p1_r = 128, p1_g = 0, p1_b = 0;
	int p2_r = 0, p2_g = 128, p2_b = 0;
	int p3_r = 0, p3_g = 0, p3_b = 128;
	int kd = 5, ka = 1;
	int ks = 100;
	double gamma = 8;
	// plane normals -- normalized
	Vect pn[] = { Vect(0,1,0), Vect(1,0,0), Vect(0,0,1) };
	double p_D[] = { 0, 0, 0 };

	// Spheres
	double s_r[] = { 0.4, 0.4, 0.4 };
	Vect sc[] = { Vect(3.0, s_r[0] + 0.01, 1.5), Vect(2.5, s_r[1] + 0.01, 2.5), Vect(1.5, s_r[2] + 0.01, 3.0) };

	// Lights
	Vect light1_pos = Vect(3, s_r[1] + 1, 3);

	Vect normal;

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
			double iamnt = (i + 0.5) / width;
			double jamnt = ((height - j) + 0.5) / height;
			Vect ray_direction = (cam_dir + cam_right*(jamnt - 0.5) +
				cam_down*(iamnt - 0.5)).normalize();
			double object_ts[OBJS_NUM];
			double object_denoms[OBJS_NUM];
			two_doubles t_d;
			// planes intersection
			for (int pl = 0; pl < 3; pl++)
			{
				intersect_plane(eye_pos, ray_direction, pn[pl], p_D[pl], t_d);
				object_ts[pl] = t_d.t;
				object_denoms[pl] = t_d.denom;
			}

			// Sphere intersections
			for (int sp = 0; sp < 3; sp++)
			{
				intersect_sphere(eye_pos, ray_direction, sc[sp], s_r[sp], t_d);
				object_denoms[sp + 3] = t_d.denom;
				object_ts[sp + 3] = t_d.t;
			}

			double denom = 0;
			double t = 0;
			int p_r, p_g, p_b;

			bool intersected_any = false;
			for (int object_t_idx = 0; object_t_idx < OBJS_NUM; object_t_idx++) {
				bool intersected = true;
				double object_t = object_ts[object_t_idx];
				if (object_t < 0 || object_t >= 1e7) {
					intersected = false;
					continue;
				}
				for (int other_object_t_idx = 0; other_object_t_idx < OBJS_NUM; other_object_t_idx++) {
					double other_object_t = object_ts[other_object_t_idx];
					if (other_object_t > 0 && other_object_t < 1e7 && other_object_t < object_t) {
						intersected = false;
						break;
					}
				}
				if (intersected) {
					intersected_any = true;
					denom = object_denoms[object_t_idx];
					t = object_ts[object_t_idx];
					if (object_t_idx < 3)
						normal = pn[object_t_idx];
					else
						normal = (eye_pos + ray_direction*t - sc[object_t_idx-3]).normalize();
					Vect intersect = eye_pos + ray_direction*t;
					double intersect_1, intersect_2;
					t_d = max(intersect.getVX(), intersect.getVY() / height * f * 2, intersect.getVZ());
					intersect_1 = t_d.denom;
					intersect_2 = t_d.t;
					if (object_t_idx < 3)
						// plane stripes
						p_r = p_g = p_b = 255 * (abs(cos(freq / f*intersect_1 * 2 * M_PI)) > thickness 
							&& abs(cos(freq / f*intersect_2 * 2 * M_PI)) > thickness);
					else {
						// spheres
						p_r = 255 * (object_t_idx == 3);
						p_g = 255 * (object_t_idx == 4);
						p_b = 255 * (object_t_idx == 5);
					}
					break;
				}
			}

			if (!intersected_any) {
				continue;
			}
			Vect intersect = eye_pos + ray_direction*t;
			Vect intersect2light = light1_pos - intersect;
			double intersect_dist = intersect2light.magnitude();
			intersect2light = intersect2light.normalize();
			bool shade = false;
			for(int pl=0; pl<3; pl++)
				if (!shade && intersect_plane(intersect, intersect2light, pn[pl], p_D[pl], t_d) && t_d.t < intersect_dist)
				{
					shade = true;
				}
			for (int sp = 0; sp<3; sp++)
				if (!shade && intersect_sphere(intersect, intersect2light, sc[sp], s_r[sp], t_d) && t_d.t < intersect_dist)
				{
					shade = true;
				}
			Vect halfway = (ray_direction*(-1) + intersect2light).normalize();
			double specular = halfway.dotProduct(normal);
			if (shade)
			{
				//cout << "shade\n";
				arr[i][j][0] = ka 
					+ kd*max(0.05, normal.dotProduct(intersect2light))*p_r / 10 / intersect_dist
					+ ks*pow(specular, gamma);
				arr[i][j][1] = ka 
					+ kd*max(0.05, normal.dotProduct(intersect2light))*p_g / 10 / intersect_dist
					+ ks*pow(specular, gamma);
				arr[i][j][2] = ka 
					+ kd*max(0.05, normal.dotProduct(intersect2light))*p_b / 10 / intersect_dist
					+ ks*pow(specular, gamma);
			}
			else
			{
				arr[i][j][0] = ka
					+ kd*max(0.05, normal.dotProduct(intersect2light))*p_r / intersect_dist / intersect_dist
					+ ks*pow(specular, gamma);
				arr[i][j][1] = ka 
					+ kd*max(0.05, normal.dotProduct(intersect2light))*p_g / intersect_dist / intersect_dist
					+ ks*pow(specular, gamma);
				arr[i][j][2] = ka 
					+ kd*max(0.05, normal.dotProduct(intersect2light))*p_b / intersect_dist / intersect_dist
					+ ks*pow(specular, gamma);
			}
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


bool intersect_plane(Vect pos, Vect ray_direction, Vect plane_normal, double D, two_doubles& td)
{
	// plane intersection
	double t;
	double denom = plane_normal.dotProduct(ray_direction);
	if (abs(denom) - 1e-5 <= 0)
		t = 1e7;
	t = -(D + plane_normal.dotProduct(pos)) / denom;
	td.t = t;
	td.denom = denom;
	return (t > 1e-5 && t < 1e7);
}

bool intersect_sphere(Vect pos, Vect ray_direction, Vect center, double radius, two_doubles& td)
{
	// Sphere intersection
	double C = ((pos - center).squaredMagnitude()) - radius*radius;
	double B = ray_direction.dotProduct(pos - center) * 2;
	double A = ray_direction.squaredMagnitude();
	double D = B*B - 4 * A * C;
	double t, t_1, t_2;
	if (D < 0) {
		// No intersection
		t = 1e7;
	}
	else if (abs(D) < 1e-5) {
		// One intersection
		t = (-1 * B) / (2 * A);
	}
	else if (D > 0) {
		// Two intersections
		t_1 = (-1 * B + sqrt(D)) / (2 * A);
		t_2 = (-1 * B - sqrt(D)) / (2 * A);
		if (t_1 < t_2)
			t = t_1;
		else
			t = t_2;
	}
	if (t > 0 && t < 1e7) {
		Vect sphere_intersect = pos + ray_direction*t;
		Vect sphere_normal = (sphere_intersect - center).normalize();
		double denom = sphere_normal.dotProduct(ray_direction);
		td.denom = denom;
	}
	td.t = t;
	return (t > 1e-5 && t < 1e7);
}