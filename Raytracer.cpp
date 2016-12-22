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
bool intersect_sphere_from_inside(Vect pos, Vect ray_direction, Vect center, double radius, two_doubles& td);

bool is_shadow(Vect start_p, Vect ray_direction, double intersect_dist, int src_obj_idx);

three_doubles compute_shade(Vect start_p, Vect ray_direction, int depth);

int width = 256 * 2, height = 256 * 2;
double freq = 200;
double thickness = 0.2;
double f = 256;
double proj_x = 3;

int global_depth = 3;
double epsilone = 0.01;

Vect O(0, 0, 0);
Vect X(1, 0, 0);
Vect Y(0, 1, 0);
Vect Z(0, 0, 1);

Vect eye_pos = Vect(5+3-3.25, 3+0.5-2.5, 5+3-3.25);
Vect look_at = Vect(0, 1, 0);
// Camera three axes are cam_dir, cam_right, cam_down
Vect cam_dir = (look_at - eye_pos).normalize();
Vect cam_right = Y.crossProduct(cam_dir).normalize();
Vect cam_down = (cam_dir.crossProduct(cam_right).normalize())*(-1);

int p1_r = 128, p1_g = 0, p1_b = 0;
int p2_r = 0, p2_g = 128, p2_b = 0;
int p3_r = 0, p3_g = 0, p3_b = 128;
int kd=3, ka=10;
int ks[] = {0,0,0,0,110,110};
double spec_gamma = 8;

float kref[] = {0.0, 0.0, 0.0, 0.0, 1.0, 0.0};
float ktrans[] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
float glass_refraction_index = 1.6, air_refraction_index = 1.0;

float gaussian_filter[] = { 0.00390625,  0.015625  ,  0.0234375 ,  0.015625  ,  0.00390625,
                            0.015625  ,  0.0625    ,  0.09375   ,  0.0625    ,  0.015625  ,
                            0.0234375 ,  0.09375   ,  0.140625  ,  0.09375   ,  0.0234375 ,
                            0.015625  ,  0.0625    ,  0.09375   ,  0.0625    ,  0.015625  ,
                            0.00390625,  0.015625  ,  0.0234375 ,  0.015625  ,  0.00390625};//{1.0/16, 1.0/8, 1.0/16, 1.0/8, 1.0/4, 1.0/8, 1.0/16, 1.0/8, 1.0/16};

// plane normals -- normalized
Vect pn[] = { Vect(0,1,0), Vect(1,0,0), Vect(0,0,1) };
double p_D[] = { 0, 0, 0 };

// Spheres
double s_r[] = { 0.4, 0.4, 0.4 };
Vect sc[] = { Vect(3.0, s_r[0] + 0.01, 1.5), Vect(2.5, s_r[1] + 0.01, 2.5), Vect(1.5, s_r[2] + 0.01, 3.0) };

// Lights
Vect light1_pos = Vect(3.5, s_r[1] + 3, 1.5); //Vect(3.5, s_r[1] + 3, 1.5)
Vect light2_pos = Vect(0.5, 0.5, 1.5); //Vect(3.5, s_r[1] + 3, 1.5)
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
            arr[i][j][0] = 0;
            arr[i][j][1] = 0;
            arr[i][j][2] = 0;

//            double x = i + 0 * 0.45;
//            double y = j + 0 * 0.45;

//            double iamnt = (x + 0.5) / width;
//            double jamnt = ((height - y) + 0.5) / height;
//            Vect ray_direction = (cam_dir + cam_right*(jamnt - 0.5) +
//                cam_down*(iamnt - 0.5)).normalize();

//            three_doubles color = compute_shade(eye_pos, ray_direction, 0);
//            arr[i][j][0] += (int)(color.r);// * gaussian_filter[index]);
//            arr[i][j][1] += (int)(color.g);// * gaussian_filter[index]);
//            arr[i][j][2] += (int)(color.b);// * gaussian_filter[index]);

            int index = 0;
            for(int horizontal = -2; horizontal <= 2; horizontal++)
                for(int vertical = -2; vertical <= 2; vertical++)
                {
                    double x = i + float((rand()%1000 - 500)/1000.0);
                    double y = j + float((rand()%1000 - 500)/1000.0);
			        
                    double iamnt = (x + 0.5) / width;
                    double jamnt = ((height - y) + 0.5) / height;
                    Vect ray_direction = (cam_dir + cam_right*(jamnt - 0.5) +
                        cam_down*(iamnt - 0.5)).normalize();
				
                    three_doubles color = compute_shade(eye_pos, ray_direction, 0);
                    arr[i][j][0] += (int)(color.r * gaussian_filter[index]);
                    arr[i][j][1] += (int)(color.g * gaussian_filter[index]);
                    arr[i][j][2] += (int)(color.b * gaussian_filter[index]);
			        
                    index++;
                }
		}
	writePPM("test.ppm", width, height, arr);
	return 0;
}

three_doubles compute_shade(Vect start_p, Vect ray_direction, int depth)
{
    three_doubles color;
    color.r = 0;
    color.g = 0;
    color.b = 0;
    if(depth == global_depth)
        return color;

    //shading
	double object_ts[OBJS_NUM];
	double object_denoms[OBJS_NUM];
	two_doubles t_d;
	// planes intersection
	for (int pl = 0; pl < 3; pl++)
	{
		intersect_plane(start_p, ray_direction, pn[pl], p_D[pl], t_d);
		object_ts[pl] = t_d.t;
		object_denoms[pl] = t_d.denom;
	}

	// Sphere intersections
	for (int sp = 0; sp < 3; sp++)
	{
		intersect_sphere(start_p, ray_direction, sc[sp], s_r[sp], t_d);
		object_denoms[sp + 3] = t_d.denom;
		object_ts[sp + 3] = t_d.t;
	}

	double denom = 0;

    double min_t = 10000000;
    int min_ind = -1;
    
    Vect normal;

    // Check shortest primary intersection
	for (int object_t_idx = 0; object_t_idx < OBJS_NUM; object_t_idx++) {
		double object_t = object_ts[object_t_idx];
		if (object_t > epsilone && object_t < 1e7 && object_t < min_t) {
			min_t = object_t;
			min_ind = object_t_idx;
		}
	}
	if(min_ind == -1)
	{
	    return color;
	}
	
	denom = object_denoms[min_ind];
	//t = object_ts[min_ind];
    if (min_ind < 3) {
		normal = pn[min_ind];
        // Gouraud normals modification
        for(int gouraud_i = 0; gouraud_i < 3; gouraud_i++){
            if( gouraud_i != min_ind ){
                double dist_diff = abs(object_ts[min_ind] - object_ts[gouraud_i]);
                if( dist_diff < 1 ){
                    normal = (pn[min_ind] + pn[gouraud_i]*(1 - 1*dist_diff)).normalize();
                    denom = normal.dotProduct(ray_direction);
                }
            }
        }
    }
	else
		normal = (start_p + ray_direction * min_t - sc[min_ind-3]).normalize();
	Vect intersect = start_p + ray_direction * min_t;
	double intersect_1, intersect_2;
	t_d = max(intersect.getVX(), intersect.getVY() / height * f * 2, intersect.getVZ());
	intersect_1 = t_d.denom;
	intersect_2 = t_d.t;
	if (min_ind < 3)
		// plane stripes
		color.r = color.g = color.b = 255 * (abs(cos(freq / f * intersect_1 * 2 * M_PI)) > thickness 
			&& abs(cos(freq / f * intersect_2 * 2 * M_PI)) > thickness);
	else {
		// spheres
        color.r = 255;
        color.g = 255; //* (min_ind == 4);
        color.b =  255;
		//color.r = color.g = color.b = 255;
	}

    Vect intersect2light1 = light1_pos - intersect;
    double intersect_dist1 = intersect2light1.magnitude();
    intersect2light1 = intersect2light1.normalize();
    Vect intersect2light2 = light2_pos - intersect;
    double intersect_dist2 = intersect2light2.magnitude();
    intersect2light2 = intersect2light2.normalize();

    bool shadow1 = is_shadow(intersect, intersect2light1, intersect_dist1, min_ind);
    bool shadow2 = is_shadow(intersect, intersect2light2, intersect_dist2, min_ind);

    Vect halfway1 = (ray_direction*(-1) + intersect2light1).normalize();
    double specular1 = halfway1.dotProduct(normal);
    Vect halfway2 = (ray_direction*(-1) + intersect2light2).normalize();
    double specular2 = halfway2.dotProduct(normal);
	
	three_doubles reflected_color;
	reflected_color.r = 0;
    reflected_color.g = 0;
    reflected_color.b = 0;

    float R = 1.0;
    
	if(kref[min_ind] > 0 || ktrans[min_ind] > 0)
	{
        Vect ref_ray_direction = ray_direction - normal * (ray_direction.dotProduct(normal)) * 2;
        ref_ray_direction = ref_ray_direction.normalize();
	    reflected_color = compute_shade(intersect, ref_ray_direction, depth + 1);
	    /*reflected_color.r /= 5;
        reflected_color.g /= 5;
        reflected_color.b /= 5;*/
    }
    
    three_doubles refracted_color;
	refracted_color.r = 0;
    refracted_color.g = 0;
    refracted_color.b = 0;
    
    if(ktrans[min_ind] > 0)
    {
        float temp1 = ray_direction.dotProduct(normal);
        float t1_b = 1 - (air_refraction_index * air_refraction_index * (1 - temp1 * temp1) / (glass_refraction_index * glass_refraction_index));
        if(t1_b >= 0)
        {
            Vect t1_ray_direction = (ray_direction - normal * temp1) * air_refraction_index / glass_refraction_index;
            t1_ray_direction = t1_ray_direction - normal * sqrt(t1_b);
            t1_ray_direction.normalize();
            R = (glass_refraction_index - 1) / (glass_refraction_index + 1);
            R *= R;
            bool existed  = intersect_sphere_from_inside(intersect, t1_ray_direction, sc[min_ind-3], s_r[min_ind-3], t_d);
            //cout<<t_d.t<<endl;
            Vect intersect_out = intersect + t1_ray_direction * t_d.t;
            Vect normal_out = (intersect_out - sc[min_ind-3]).normalize() * -1;
            float temp2 = t1_ray_direction.dotProduct(normal_out);
            float t2_b = 1 - (glass_refraction_index * glass_refraction_index * (1 - temp2 * temp2) / (air_refraction_index * air_refraction_index));
            if(t2_b >= 0)
            {
               Vect t2_ray_direction = (t1_ray_direction - normal_out * temp2) * glass_refraction_index / air_refraction_index;
               t2_ray_direction = t2_ray_direction - normal_out * sqrt(t2_b);
               t2_ray_direction.normalize();         
               refracted_color = compute_shade(intersect_out, t2_ray_direction, depth + 1);
            }
        }
    }
    double shadow1_atten = 1;
    double shadow2_atten = 1;
    if (shadow1)
        shadow1_atten = 5;
    if (shadow2)
        shadow2_atten = 5;

    color.r = ka
        + kd*max(0.05, normal.dotProduct(intersect2light1))*color.r / intersect_dist1/ intersect_dist1/ shadow1_atten
        + ks[min_ind]*pow(specular1, spec_gamma)
        + kd*max(0.05, normal.dotProduct(intersect2light2))*color.r / intersect_dist2/ intersect_dist2/ shadow2_atten
        + ks[min_ind]*pow(specular2, spec_gamma)
        + kref[min_ind] * reflected_color.r
        + ktrans[min_ind] * ((1 - R) *refracted_color.r + R * reflected_color.r);

    color.g = ka
        + kd*max(0.05, normal.dotProduct(intersect2light1))*color.g / intersect_dist1 / intersect_dist1/ shadow1_atten
        + ks[min_ind]*pow(specular1, spec_gamma)
        + kd*max(0.05, normal.dotProduct(intersect2light2))*color.g / intersect_dist2/ intersect_dist2/ shadow2_atten
        + ks[min_ind]*pow(specular2, spec_gamma)
        + kref[min_ind] * reflected_color.g
        + ktrans[min_ind] * ((1 - R) *refracted_color.g + R * reflected_color.g);

    color.b = ka
        + kd*max(0.05, normal.dotProduct(intersect2light1))*color.b / intersect_dist1 / intersect_dist1/ shadow1_atten
        + ks[min_ind]*pow(specular1, spec_gamma)
        + kd*max(0.05, normal.dotProduct(intersect2light2))*color.b / intersect_dist2/ intersect_dist2/ shadow2_atten
        + ks[min_ind]*pow(specular2, spec_gamma)
        + kref[min_ind] * reflected_color.b
        + ktrans[min_ind] * ((1 - R) *refracted_color.b + R * reflected_color.b);

    //color clipping    
    color.r = color.r < 0? 0: color.r;
    color.g = color.g < 0? 0: color.g;
    color.b = color.b < 0? 0: color.b;  

    color.r = color.r > 255? 255: color.r;
    color.g = color.g > 255? 255: color.g;
    color.b = color.b > 255? 255: color.b;  
    
    return color;
}

bool is_shadow(Vect start_p, Vect ray_direction, double intersect_dist, int src_obj_idx)
{
	two_doubles t_d;
	bool shade = false;
	for(int pl=0; pl<3; pl++)
        if (!shade && intersect_plane(start_p, ray_direction, pn[pl], p_D[pl], t_d) && (src_obj_idx != pl) && t_d.t < intersect_dist)
			shade = true;
	for (int sp = 0; sp<3; sp++)
        if (!shade && intersect_sphere(start_p, ray_direction, sc[sp], s_r[sp], t_d) && (src_obj_idx - 3 != sp) && t_d.t < intersect_dist)
			shade = true;
    return shade;
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

bool intersect_sphere_from_inside(Vect pos, Vect ray_direction, Vect center, double radius, two_doubles& td)
{
	// Sphere intersection
	double C = ((pos - center).squaredMagnitude()) - radius*radius;
	double B = ray_direction.dotProduct(pos - center) * 2;
	double A = ray_direction.squaredMagnitude();
	double D = B*B - 4 * A * C;
	double t, t_1, t_2;
	//cout<<D<<endl;
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
		if (t_1 > t_2)
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
