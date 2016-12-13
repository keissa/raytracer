//#define _USE_MATH_DEFINES
//#include <cmath>
//#include <fstream>
//#include <iostream>

//#define OBJS_NUM 1

//using namespace std;

//void writePPM(const char * imagepath, int width, int height, int*** img, int color = 255);

//int main()
//{
//    int width = 256*2, height = 256*2;
//    double freq = 3;
//    double thickness = 0.2;
////    double f = 15;
//    double f = 256;
//    double proj_x = 3;
////    double e_x = 25, e_y = 0, e_z = 25;
//    double e_x = 0, e_y = 0, e_z = 512;
////    int p1_r = 128, p1_g = 0, p1_b = 0;
////    int p2_r = 0, p2_g = 128, p2_b = 0;
////    int p3_r = 0, p3_g = 0, p3_b = 128;
//    int kd = 40, ka = 10;
//    // plane normal -- normalized
//    double angle = 45;
//    double n_x = cos(M_PI * angle / 180), n_y = 0, n_z = sin(M_PI * angle / 180), D = -f * sin(M_PI * angle / 180);

//    double sc1_x = 0, sc1_y = 0, sc1_z = 0, s1_r = 100;

//    int*** arr = new int**[width];
//    for (int i = 0; i < width; i++)
//    {
//        arr[i] = new int*[height];
//        for (int j = 0; j < height; j++)
//            arr[i][j] = new int[3]();
//    }
//    for (int i = 0; i < width; i++)
//        for (int j = 0; j < height; j++)
//        {
////            double x = (((j - 0)*(f - proj_x - proj_x)) / (height - 1)) + proj_x;

//            double x = j - width/ 2.0;
////            double y = height*0.75-i;
////            double y = 4.5-i*9/width;
//            double y = -i + height/ 2.0;
//            double dd_x = x - e_x, dd_y = y - e_y;
////            double z = -1*(x*n_x + y*n_y + D)/n_z;
//            double z = f;
//            double dd_z = z - e_z;
//            double norm = sqrt(dd_x*dd_x + dd_y*dd_y + dd_z*dd_z);
//            double d_x = dd_x / norm;
//            double d_y = dd_y / norm;
//            double d_z = dd_z / norm;

//            double object_ts[OBJS_NUM];
//            double object_denoms[OBJS_NUM];

//            int p_r, p_g, p_b;
//            // Sphere intersection
//            // Right Sphere
//            //
//            double C = ( (e_x-sc1_x)*(e_x-sc1_x) + (e_y-sc1_y)*(e_y-sc1_y) + (e_z-sc1_z)*(e_z-sc1_z) ) - s1_r*s1_r;
//            double B = 2 * (d_x*(e_x-sc1_x) + d_y*(e_y-sc1_y) + d_z*(e_z-sc1_z));
//            double A = d_x*d_x + d_y*d_y + d_z*d_z;
//            double D = B*B - 4 * A * C;
//            double t4, t4_1, t4_2;
//            if ( D < 0 ){
//                // No intersection
//                t4 = 1e7;
//            } else if ( abs(D) < 1e-2 ) {
//                // One intersection
//                cout << "1" << endl;
//                t4 = (-1*B) / (2*A);
//            } else if ( D > 0 ) {
//                // Two intersections

//                t4_1 = (-1*B + sqrt(D)) / (2*A);
//                t4_2 = (-1*B - sqrt(D)) / (2*A);
//                if( t4_1 < t4_2 )
//                    t4 = t4_1;
//                else
//                    t4 = t4_2;
////                cout << "2 " << t4 << " " << t4_1 << " " << t4_2 << endl;
//            }
//            if( t4 > 0 && t4 < 1e7 ){

//                double s_intersect_x = e_x + t4*d_x;
//                double s_intersect_y = e_y + t4*d_y;
//                double s_intersect_z = e_z + t4*d_z;

//                double sn1_x = s_intersect_x - sc1_x, sn1_y = s_intersect_y - sc1_y, sn1_z = s_intersect_z - sc1_z;
//                double norm = sqrt(sn1_x*sn1_x + sn1_y*sn1_y + sn1_z*sn1_z);
//                sn1_x /= norm;
//                sn1_y /= norm;
//                sn1_z /= norm;
//                double denom4 = (sn1_x*d_x + sn1_y*d_y + sn1_z*d_z);
//                object_denoms[0] = denom4;
////                cout << "coloring red " << denom4 << " " << endl;
//                p_r = 255;
//                p_g = 0;
//                p_b = 0;
//            } else {
//                p_r = 0;
//                p_g = 0;
//                p_b = 0;
//            }
//            object_ts[0] = t4;

////            cout << "asdasd " << ka + kd*abs(object_denoms[0])*p_r << endl;
//            arr[i][j][0] = ka + kd*abs(object_denoms[0])*p_r;
//            arr[i][j][1] = ka + kd*abs(object_denoms[0])*p_g;
//            arr[i][j][2] = ka + kd*abs(object_denoms[0])*p_b;

//            // clip values
//            arr[i][j][0] = arr[i][j][0] > 255 ? 255 : arr[i][j][0];
//            arr[i][j][1] = arr[i][j][1] > 255 ? 255 : arr[i][j][1];
//            arr[i][j][2] = arr[i][j][2] > 255 ? 255 : arr[i][j][2];
//        }
//    writePPM("test.ppm", width, height, arr);
//    return 0;
//}

//void writePPM(const char * imagepath, int width, int height, int*** img, int color)
//{
//    std::ofstream f(imagepath);
//    // header
//    f << "P3" << "\n";
//    // meta
//    f << width << " " << height << " " << color << "\n";
//    // data
//    int r = 255, g = 0, b = 0;
//    for (int i = 0; i < width; i++)
//        for (int j = 0; j < height; j++)
//            f << img[i][j][0] << " " << img[i][j][1] << " " << img[i][j][2] << "\n";

//    f.close();
//}
