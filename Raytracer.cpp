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

int main()
{
	int p1_r = 128, p1_g = 0, p1_b = 0;
	int p2_r = 0, p2_g = 128, p2_b = 0;
	int p3_r = 0, p3_g = 0, p3_b = 128;
	int kd = 5, ka = 1;
    int ks = 200;
	double gamma = 8;
	// plane normal -- normalized
	Vect pn1 = Vect(0, 1, 0);
	double p1_D = 0;
	Vect pn2 = Vect(1, 0, 0);
	double p2_D = 0;
	Vect pn3 = Vect(0, 0, 1);
	double p3_D = 0;

	// Right Sphere
	double s1_r = 0.4;
	Vect sc1 = Vect(3.0, s1_r + 0.01, 1.5);
	// Center Sphere
	double s2_r = 0.4;
	Vect sc2 = Vect(2.5, s2_r + 0.01, 2.5);
	// Left Sphere
	double s3_r = 0.4;
	Vect sc3 = Vect(1.5, s3_r + 0.01, 3.0);

	// Lights
	Vect light1_pos = Vect(3, s2_r + 1, 3);

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
			intersect_plane(eye_pos, ray_direction, pn1, p1_D, t_d);
			object_ts[0] = t_d.t;
			object_denoms[0] = t_d.denom;
			intersect_plane(eye_pos, ray_direction, pn2, p2_D, t_d);
			object_ts[1] = t_d.t;
			object_denoms[1] = t_d.denom;
			intersect_plane(eye_pos, ray_direction, pn3, p3_D, t_d);
			object_ts[2] = t_d.t;
			object_denoms[2] = t_d.denom;

			// Sphere intersection
			// Right Sphere
			intersect_sphere(eye_pos, ray_direction, sc1, s1_r, t_d);
			object_denoms[3] = t_d.denom;
			object_ts[3] = t_d.t;
			// Center Sphere
			intersect_sphere(eye_pos, ray_direction, sc2, s2_r, t_d);
			object_denoms[4] = t_d.denom;
			object_ts[4] = t_d.t;
			// Left Sphere
			intersect_sphere(eye_pos, ray_direction, sc3, s3_r, t_d);
			object_denoms[5] = t_d.denom;
			object_ts[5] = t_d.t;

			double denom = 0;
			double t = 0;
			int p_r, p_g, p_b;

			bool intersected_any = false;
            // Check shortest primary intersection
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
                    Vect intersect = eye_pos + ray_direction*t;
					if (object_t_idx == 0)
						normal = pn1;
					else if (object_t_idx == 1)
						normal = pn2;
					else if (object_t_idx == 2)
						normal = pn3;
					else if (object_t_idx == 3)
                        normal = (intersect - sc1).normalize();
					else if (object_t_idx == 4)
                        normal = (intersect - sc2).normalize();
					else if (object_t_idx == 5)
                        normal = (intersect - sc3).normalize();
					double intersect_1, intersect_2;
					if (object_t_idx == 0) {
						intersect_1 = intersect.getVX();
						intersect_2 = intersect.getVZ();
					}
					else if (object_t_idx == 1) {
						intersect_1 = intersect.getVZ();
						intersect_2 = intersect.getVY() / height * f * 2;
					}
					else if (object_t_idx == 2) {
						intersect_1 = intersect.getVX();
						intersect_2 = intersect.getVY() / height * f * 2;
					}
					// plane stripes
					if (object_t_idx < 3) {
						if (abs(cos(freq / f*intersect_1 * 2 * M_PI)) < thickness || abs(cos(freq / f*intersect_2 * 2 * M_PI)) < thickness)
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
					}
					else if (object_t_idx == 3) {
						// Right Sphere
						p_r = 255;
						p_g = 0;
						p_b = 0;
					}
					else if (object_t_idx == 4) {
						// Center Sphere
						p_r = 0;
						p_g = 255;
						p_b = 0;
					}
					else if (object_t_idx == 5) {
						// Left Sphere
						p_r = 0;
						p_g = 0;
						p_b = 255;
					}

					break;
				}
			}

			if (!intersected_any) {
				continue;
			}
            // Check for shadows
			Vect intersect = eye_pos + ray_direction*t;
			Vect intersect2light = light1_pos - intersect;
			double intersect_dist = intersect2light.magnitude();
			intersect2light = intersect2light.normalize();
			bool shade = false;
			if (intersect_plane(intersect, intersect2light, pn1, p1_D, t_d) && t_d.t < intersect_dist)
				shade = true;
			else if (intersect_plane(intersect, intersect2light, pn2, p2_D, t_d) && t_d.t < intersect_dist)
				shade = true;
			else if (intersect_plane(intersect, intersect2light, pn3, p3_D, t_d) && t_d.t < intersect_dist)
				shade = true;
			else if (intersect_sphere(intersect, intersect2light, sc1, s1_r, t_d) && t_d.t < intersect_dist)
				shade = true;
			else if (intersect_sphere(intersect, intersect2light, sc2, s3_r, t_d) && t_d.t < intersect_dist)
				shade = true;
			else if (intersect_sphere(intersect, intersect2light, sc3, s3_r, t_d) && t_d.t < intersect_dist)
				shade = true;
            Vect halfway = (ray_direction*(-1) + intersect2light).normalize();
			double specular = halfway.dotProduct(normal);

            // Check for first secondary reflection
            // has equation: ray2 = intersect_point + t*ray2_direction

            Vect sec_ref1_direction = ray_direction - normal * (ray_direction.dotProduct(normal)) * 2; // direction of reflected
            sec_ref1_direction = sec_ref1_direction.normalize();
            Vect second_normal;
            double ref1_t;
            double ref1_a = ka, ref1_dr = 0, ref1_dg = 0, ref1_db = 0, ref1_s = 0, ref1_coeff = 0.5;

            // planes intersection
            intersect_plane(intersect, sec_ref1_direction, pn1, p1_D, t_d);
            object_ts[0] = t_d.t;
            object_denoms[0] = t_d.denom;
            intersect_plane(intersect, sec_ref1_direction, pn2, p2_D, t_d);
            object_ts[1] = t_d.t;
            object_denoms[1] = t_d.denom;
            intersect_plane(intersect, sec_ref1_direction, pn3, p3_D, t_d);
            object_ts[2] = t_d.t;
            object_denoms[2] = t_d.denom;

            // Sphere intersection
            // Right Sphere
            intersect_sphere(intersect, sec_ref1_direction, sc1, s1_r, t_d);
            object_denoms[3] = t_d.denom;
            object_ts[3] = t_d.t;
            // Center Sphere
            intersect_sphere(intersect, sec_ref1_direction, sc2, s2_r, t_d);
            object_denoms[4] = t_d.denom;
            object_ts[4] = t_d.t;
            // Left Sphere
            intersect_sphere(intersect, sec_ref1_direction, sc3, s3_r, t_d);
            object_denoms[5] = t_d.denom;
            object_ts[5] = t_d.t;

//            double denom = 0;
//            double t = 0;
            int p_r_2, p_g_2, p_b_2;

            intersected_any = false;

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
                    Vect ref1_intersect = intersect + sec_ref1_direction*t;
                    Vect ref1_intersect2light = light1_pos - ref1_intersect;
                    double ref1_intersect_dist = ref1_intersect2light.magnitude();
                    ref1_intersect2light = ref1_intersect2light.normalize();
                    if (object_t_idx == 0)
                        second_normal = pn1;
                    else if (object_t_idx == 1)
                        second_normal = pn2;
                    else if (object_t_idx == 2)
                        second_normal = pn3;
                    else if (object_t_idx == 3)
                        second_normal = (ref1_intersect - sc1).normalize();
                    else if (object_t_idx == 4)
                        second_normal = (ref1_intersect - sc2).normalize();
                    else if (object_t_idx == 5)
                        second_normal = (ref1_intersect - sc3).normalize();
                    double intersect_1, intersect_2;
                    if (object_t_idx == 0) {
                        intersect_1 = ref1_intersect.getVX();
                        intersect_2 = ref1_intersect.getVZ();
                    }
                    else if (object_t_idx == 1) {
                        intersect_1 = ref1_intersect.getVZ();
                        intersect_2 = ref1_intersect.getVY() / height * f * 2;
                    }
                    else if (object_t_idx == 2) {
                        intersect_1 = ref1_intersect.getVX();
                        intersect_2 = ref1_intersect.getVY() / height * f * 2;
                    }
                    // plane stripes
                    if (object_t_idx < 3) {
                        if (abs(cos(freq / f*intersect_1 * 2 * M_PI)) < thickness || abs(cos(freq / f*intersect_2 * 2 * M_PI)) < thickness)
                        {
                            p_r_2 = 0;
                            p_g_2 = 0;
                            p_b_2 = 0;
                        }
                        else
                        {
                            p_r_2 = 255;
                            p_g_2 = 255;
                            p_b_2 = 255;
                        }
                    }
                    else if (object_t_idx == 3) {
                        // Right Sphere
                        p_r_2 = 255;
                        p_g_2 = 0;
                        p_b_2 = 0;
                    }
                    else if (object_t_idx == 4) {
                        // Center Sphere
                        p_r_2 = 0;
                        p_g_2 = 255;
                        p_b_2 = 0;
                    }
                    else if (object_t_idx == 5) {
                        // Left Sphere
                        p_r_2 = 0;
                        p_g_2 = 0;
                        p_b_2 = 255;
                    }
                    ref1_dr = kd*max(0.05, second_normal.dotProduct(ref1_intersect2light))*p_r_2 / 10 / ref1_intersect_dist;

                    ref1_dg = kd*max(0.05, second_normal.dotProduct(ref1_intersect2light))*p_g_2 / 10 / ref1_intersect_dist;

                    ref1_db = kd*max(0.05, second_normal.dotProduct(ref1_intersect2light))*p_b_2 / 10 / ref1_intersect_dist;
                    Vect halfway2 = (sec_ref1_direction*(-1) + ref1_intersect2light).normalize();
                    double specular2 = halfway2.dotProduct(second_normal);
                    ref1_s = ks*pow(specular2, gamma);

                    break;
                }

                if (!intersected_any) {
                    ref1_coeff = 0;
                }
            }

			if (shade)
			{
				//cout << "shade\n";
				arr[i][j][0] = ka 
					+ kd*max(0.05, normal.dotProduct(intersect2light))*p_r / 10 / intersect_dist
                    + ks*pow(specular, gamma)
                    + ref1_coeff*(ref1_a + ref1_dr + ref1_s);
				arr[i][j][1] = ka 
					+ kd*max(0.05, normal.dotProduct(intersect2light))*p_g / 10 / intersect_dist
                    + ks*pow(specular, gamma)
                    + ref1_coeff*(ref1_a + ref1_dg + ref1_s);
				arr[i][j][2] = ka 
					+ kd*max(0.05, normal.dotProduct(intersect2light))*p_b / 10 / intersect_dist
                    + ks*pow(specular, gamma)
                    + ref1_coeff*(ref1_a + ref1_db + ref1_s);
			}
			else
			{
				arr[i][j][0] = ka
					+ kd*max(0.05, normal.dotProduct(intersect2light))*p_r / intersect_dist / intersect_dist
                    + ks*pow(specular, gamma)
                    + ref1_coeff*(ref1_a + ref1_dr + ref1_s);
				arr[i][j][1] = ka 
					+ kd*max(0.05, normal.dotProduct(intersect2light))*p_g / intersect_dist / intersect_dist
                    + ks*pow(specular, gamma)
                    + ref1_coeff*(ref1_a + ref1_dg + ref1_s);
				arr[i][j][2] = ka 
					+ kd*max(0.05, normal.dotProduct(intersect2light))*p_b / intersect_dist / intersect_dist
                    + ks*pow(specular, gamma)
                    + ref1_coeff*(ref1_a + ref1_db + ref1_s);
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
