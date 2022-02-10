#include <bits/stdc++.h>
#include <windows.h>
#include <glut.h>
#define pi (2*acos(0.0))
using namespace std;
const double eps=0.0000001;
const double max_value=9999999999999;
extern int recursion_level;

struct vector3D
{
    double x,y,z;

    void set_point(double i, double  j, double k)
    {
        x = i;
        y = j;
        z = k;
    }
    double magnitude() const
    {
        return sqrt(x*x + y*y + z*z);
    }
    void normalize()
    {
        double temp = magnitude();
        if(temp<eps)
        {
            x = 0;
            y = 0;
            z = 0;
            return;
        }
        x=x/temp, y=y/temp, z=z/temp;
    }
    struct vector3D operator + (struct vector3D const &point) const
    {
        struct vector3D res{};
        res.x = x + point.x;
        res.y = y + point.y;
        res.z = z + point.z;
        return res;
    }
    struct vector3D operator - (struct vector3D const &point) const
    {
        struct vector3D res{};
        res.x = x - point.x;
        res.y = y - point.y;
        res.z = z - point.z;
        return res;
    }
    struct vector3D operator * (double d) const
    {
        struct vector3D res{};
        res.x = x * d;
        res.y = y * d;
        res.z = z * d;
        return res;
    }
    struct vector3D operator / (double d) const
    {
        struct vector3D res{};
        res.x = x / d;
        res.y = y / d;
        res.z = z / d;
        return res;
    }
    void print_point() const
    {
        cout<<"x: "<<x<<" y: "<<y<<" z: "<<z<<endl;
    }
};

double dot_product(struct vector3D a, struct vector3D b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
struct vector3D cross_product(struct vector3D a, struct vector3D b)
{
    struct vector3D res{};

    res.set_point(a.y * b.z - a.z * b.y,a.z * b.x - a.x * b.z,a.x * b.y - a.y * b.x);
    return res;
}

class Ray{
public:
    vector3D start{};
    vector3D dir{};

    Ray(vector3D s, vector3D d){
        start=s;
        dir=d;
        dir.normalize();
    }
    vector3D normalized_dir()
    {
        double temp = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
        dir.set_point(dir.x/temp, dir.y/temp, dir.z/temp);
        return dir;
    }
};

class Light{
public:
    vector3D light_pos{};
    double color[3]{};
    Light(vector3D pos, double r, double g, double b)
    {
        light_pos=pos;
        color[0]=r;
        color[1]=g;
        color[2]=b;
    }
    void draw()
    {
        glColor3f(color[0],color[1],color[2]);
        glBegin(GL_QUADS);{
            glVertex3f( light_pos.x, light_pos.y-2,light_pos.z-2);
            glVertex3f( light_pos.x,light_pos.y-2,light_pos.z+2);
            glVertex3f(light_pos.x,light_pos.y+2,light_pos.z+2);
            glVertex3f(light_pos.x, light_pos.y+2,light_pos.z-2);
        }glEnd();

    }
};

extern vector<Light*> lights;
class Object;
extern vector<Object*> objects;

class Object{
public:
    vector3D ref_point{};
    double height, width, length;
    double color[3]{};
    double coEfficients[4]{};
    int shine;
    string obj_type;
    double eqn_val[10]{};
    vector3D tri_points[3]{};
    double eta;

    Object()
    {
        eta=0.0;
        obj_type="";
        ref_point.set_point(0.0,0.0,0.0);
        height=0.0;
        width=0.0;
        length=0.0;
        color[0]=0.0;
        color[1]=0.0;
        color[2]=0.0;
        coEfficients[0]=0.0;
        coEfficients[1]=0.0;
        coEfficients[2]=0.0;
        coEfficients[3]=0.0;
        shine=0.0;
    }
    virtual void draw()
    {

    }
    void setColor(double r, double g, double b)
    {
        color[0]=r;
        color[1]=g;
        color[2]=b;
    }
    void setShine(int s)
    {
        shine=s;
    }
    void setCoefficient(double a, double b, double c, double d)
    {
        coEfficients[0]=a;
        coEfficients[1]=b;
        coEfficients[2]=c;
        coEfficients[3]=d;
    }
    void set_eqn_val(vector<double> arr)
    {
        for (int i = 0; i < 10; ++i) {
            eqn_val[i]=arr[i];
        }
    }
    virtual double intersect(Ray *r, double *newColor, int level)
    {
        return -1.0;
    }
    void lighting(double *newColor,const double *intersection_color, vector3D intersection_point,vector3D normal, Ray *r, int level) const
    {
        vector3D Rd = r->dir;
        double aux = dot_product(Rd*-1,normal);
        if(aux<0){
            normal = normal * -1;
        }
        newColor[0]=intersection_color[0]*coEfficients[0];
        newColor[1]=intersection_color[1]*coEfficients[0];
        newColor[2]=intersection_color[2]*coEfficients[0];

        unsigned int no_of_lights =lights.size();
        unsigned int no_of_objects =objects.size();

        for(auto lit : lights)
        {
            vector3D ray_dir = lit->light_pos - intersection_point;
            Ray ray_l(intersection_point,ray_dir);
            ray_l.start = ray_l.start + ray_dir*eps;

            double shadow=ray_dir.magnitude();
            double dummy_color[3];
            bool shadow_flag = false;
            for(int j=0;j<no_of_objects;j++)
            {
                double new_t = objects[j]->intersect(&ray_l,dummy_color,0);
                if(0.0<new_t && new_t<shadow)
                {
                    shadow_flag = true;
                    break;
                }
            }
            if(shadow_flag)
            {
                continue;
            }
            double theta, lambert_value, phong_value;
            lambert_value = max(dot_product(ray_l.dir,normal),0.0);
            theta = 2*dot_product(normal,ray_l.dir);

            vector3D new_v{0,0,0};
            new_v = new_v - Rd;

            vector3D temp_r = normal * theta;
            temp_r = temp_r - ray_l.dir;
            Ray ray_r(intersection_point,temp_r);

            phong_value = max(dot_product(ray_r.dir,new_v),0.0);

            for(int j=0;j<3;j++)
            {
                newColor[j]+=lit->color[j]*coEfficients[1]*lambert_value*intersection_color[j];
                newColor[j]+=lit->color[j]*coEfficients[2]*pow(phong_value,shine);
                if(newColor[j]<0) newColor[j]=0.0;
                else if(newColor[j]>1) newColor[j]=1.0;
            }
        }

        if(level>=recursion_level)
        {
            return;
        }
        double rlf_color[3];
        int nearest=-1;
        double t=max_value;

        vector3D ref = Rd - normal* 2* dot_product(normal,Rd);
        Ray ref_ray(intersection_point,ref);
        ref_ray.start = ref_ray.start + ref_ray.dir*eps;

        for(int k=0; k<no_of_objects; k++)
        {
            double v = objects[k]->intersect(&ref_ray,rlf_color,0);

            if(0<v && v<t){
                t = v;
                nearest = k ;
            }
        }
        if(nearest!=-1){
            objects[nearest]->intersect(&ref_ray,rlf_color,level+1);
            for (int j=0; j<3; j++) {
                newColor[j]+=rlf_color[j]*coEfficients[3];
                if(newColor[j]<0) newColor[j]=0.0;
                else if(newColor[j]>1) newColor[j]=1.0;
            }
        }

        if(abs(eta)<eps)
        {
            return;
        }
        double refract_color[3];
        nearest=-1;
        t=max_value;

        //double eta=0.3;
        vector3D refract{};
        double temp1 = dot_product(normal,Rd);
        double temp2 = 1 - eta * eta * (1 - temp1*temp1);
        if(temp2< 0)
        {
            return;
            //refract.set_point(0,0,0);
        }
        else
        {
            refract = Rd * eta - normal * (eta * temp1 + sqrt(temp2));
        }
        Ray refract_ray(intersection_point,refract);
        refract_ray.start = refract_ray.start + refract_ray.dir*eps;

        for(int k=0; k<no_of_objects; k++)
        {
            double v = objects[k]->intersect(&refract_ray,refract_color,0);

            if(0<v && v<t){
                t = v;
                nearest = k ;
            }
        }
        if(nearest!=-1){
            objects[nearest]->intersect(&refract_ray,refract_color,level+1);
            for (int j=0; j<3; j++) {
                newColor[j]+=refract_color[j];
                if(newColor[j]<0) newColor[j]=0.0;
                else if(newColor[j]>1) newColor[j]=1.0;
            }
        }
    }
};


class Quadric : public Object{
public:
    Quadric(vector3D point, double l, double w, double h)
    {
        //eta = 0.1;
        ref_point = point;
        length = l;
        width = w;
        height = h;
    }

    bool is_inside_cube(vector3D intersection_point)
    {
        return (abs(length)<eps || (ref_point.x<intersection_point.x && intersection_point.x<ref_point.x+length)) &&
               (abs(width)<eps || (ref_point.y<intersection_point.y && intersection_point.y<ref_point.y+width)) &&
               (abs(height)<eps || (ref_point.z<intersection_point.z && intersection_point.z<ref_point.z+height)) ;
    }

    double intersect(Ray *r, double *newColor, int level) override
    {
        vector3D Ro = r->start, Rd = r->dir;
        double a,b,c,d,temp,t_pos,t_neg,t_min;

        a = eqn_val[0]*Rd.x*Rd.x + eqn_val[1]*Rd.y*Rd.y + eqn_val[2]*Rd.z*Rd.z + eqn_val[3]*Rd.x*Rd.y + eqn_val[4]*Rd.x*Rd.z + eqn_val[5]*Rd.y*Rd.z;
        b = 2*eqn_val[0]*Ro.x*Rd.x + 2*eqn_val[1]*Ro.y*Rd.y + 2*eqn_val[2]*Ro.z*Rd.z + eqn_val[3]*Ro.x*Rd.y + eqn_val[3]*Rd.x*Ro.y +
            eqn_val[4]*Ro.x*Rd.z + eqn_val[4]*Ro.z*Rd.x + eqn_val[5]*Ro.z*Rd.y + eqn_val[5]*Ro.y*Rd.z + eqn_val[6]*Rd.x + eqn_val[7]*Rd.y + eqn_val[8]*Rd.z;
        c = eqn_val[0]*Ro.x*Ro.x + eqn_val[1]*Ro.y*Ro.y + eqn_val[2]*Ro.z*Ro.z + eqn_val[3]*Ro.x*Ro.y + eqn_val[4]*Ro.x*Ro.z + eqn_val[5]*Ro.z*Ro.y +
            eqn_val[6]*Ro.x + eqn_val[7]*Ro.y + eqn_val[8]*Ro.z + eqn_val[9];

        if(abs(a)<eps)
        {
            if(abs(b)<eps) return -1;

            t_min=-c/b;
            if(t_min<=0 || !is_inside_cube(Ro+Rd*t_min)) return -1;
        }
        else {
            temp = b * b - 4 * a * c;
            if (temp < 0.0) {
                t_min = -1.0;
                return t_min;
            } else {
                d = sqrt(temp);
                t_neg = (-b - d) / (2 * a);
                t_pos = (-b + d) / (2 * a);

                if (t_pos < t_neg) {
                    swap(t_pos, t_neg);
                }
                if (t_neg > 0 && is_inside_cube(Ro + Rd * t_neg)) {
                    t_min = t_neg;
                } else if (t_pos > 0 && is_inside_cube(Ro + Rd * t_pos)) {
                    t_min = t_pos;
                } else {
                    return -1.0;
                }
            }
        }
        if(level==0) return t_min;

        vector3D intersection_point{}, normal{};
        intersection_point = r->start + r->dir * t_min;
        normal.set_point(2*eqn_val[0]*intersection_point.x+eqn_val[3]*intersection_point.y+eqn_val[4]*intersection_point.z+eqn_val[6],
                         2*eqn_val[1]*intersection_point.y+eqn_val[3]*intersection_point.x+eqn_val[5]*intersection_point.z+eqn_val[7],
                         2*eqn_val[2]*intersection_point.z+eqn_val[4]*intersection_point.x+eqn_val[5]*intersection_point.y+eqn_val[8]);

        normal.normalize();

        lighting(newColor, color, intersection_point, normal, r, level) ;

        return t_min;
    }
};

class Sphere : public Object{
public:
    Sphere(vector3D center, double radius)
    {
        //eta = 0.3;
        ref_point=center;
        length=radius;
    }
    void draw() override
    {
        struct vector3D points[100][100];
        int i,j;
        double h,r;
        int stacks=90, slices=50;
        //generate points
        for(i=0;i<=stacks;i++)
        {
            h=length*sin(((double)i/(double)stacks)*(pi/2));
            r=length*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0;j<=slices;j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }
        //draw quads using generated points
        for(i=0;i<stacks;i++)
        {
            glColor3f(color[0],color[1],color[2]);
            for(j=0;j<slices;j++)
            {
                glBegin(GL_QUADS);{
                    //upper hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                    //lower hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                }glEnd();
            }
        }

    }

    double intersect(Ray *r, double *newColor, int level) override
    {
        double a,b,c,d,t_pos,t_neg,t_min,temp;
        vector3D Ro{}, Rd{};
        Ro = r->start - ref_point;
        Rd = r->dir;
        a = 1.0;
        b = 2 * dot_product(Ro,Rd);
        c = dot_product(Ro,Ro) - length*length;
        temp = b*b-4*a*c;

        if(temp<0)
        {
            t_min=-1.0;
            return t_min;
        }
        else{
            d = sqrt(temp);
            t_pos = (-b + d)/(2*a);
            t_neg = (-b - d)/(2*a);
            if(t_pos<0)
            {
                t_min=-1.0;
                return t_min;
            }
            else if(t_pos>0 && t_neg<0)
            {
                t_min=t_pos;
            }
            else{
                t_min=t_neg;
//                if(c>0){
//                    t_min = t_neg;
//                }
//                else {
//                    t_min = t_pos;
//                }
            }
        }
        if(level==0) return t_min;

        vector3D intersection_point{}, normal{};
        intersection_point = r->start + r->dir * t_min;
        normal = intersection_point - ref_point;
        normal.normalize();

        lighting(newColor, color, intersection_point, normal, r, level) ;


//        newColor[0]=color[0]*coEfficients[0];
//        newColor[1]=color[1]*coEfficients[0];
//        newColor[2]=color[2]*coEfficients[0];
//
//        unsigned int no_of_lights =lights.size();
//        unsigned int no_of_objects =objects.size();
//
//        for(auto lit : lights)
//        {
//            vector3D ray_dir = lit->light_pos - intersection_point;
//            Ray ray_l(intersection_point,ray_dir);
//            ray_l.start = ray_l.start + ray_dir*eps;
//
//            double shadow=ray_dir.magnitude();
//            double dummy_color[3];
//            bool shadow_flag = false;
//            for(int j=0;j<no_of_objects;j++)
//            {
//                double new_t = objects[j]->intersect(&ray_l,dummy_color,0);
//                if(0.0<new_t && new_t<shadow)
//                {
//                    shadow_flag = true;
//                    break;
//                }
//            }
//            if(shadow_flag)
//            {
//                continue;
//            }
//            double theta, lambert_value, phong_value;
//            lambert_value = max(dot_product(ray_l.dir,normal),0.0);
//            theta = 2*dot_product(normal,ray_l.dir);
//
//            vector3D new_v{0,0,0};
//            new_v = new_v - Rd;
//
//            vector3D temp_r = normal * theta;
//            temp_r = temp_r - ray_l.dir;
//            Ray ray_r(intersection_point,temp_r);
//
//            phong_value = max(dot_product(ray_r.dir,new_v),0.0);
//
//            for(int j=0;j<3;j++)
//            {
//                newColor[j]+=lit->color[j]*coEfficients[1]*lambert_value*color[j];
//                newColor[j]+=lit->color[j]*coEfficients[2]*pow(phong_value,shine);
//                if(newColor[j]<0) newColor[j]=0.0;
//                else if(newColor[j]>1) newColor[j]=1.0;
//            }
//        }
//
//        if(level>=recursion_level)
//        {
//            return t_min;
//        }
//        double rlf_color[3];
//        int nearest=-1;
//        double t=max_value;
//
//        vector3D ref = Rd - normal* 2* dot_product(normal,Rd);
//        Ray ref_ray(intersection_point,ref);
//        ref_ray.start = ref_ray.start + ref_ray.dir*eps;
//
//        for(int k=0; k<no_of_objects; k++)
//        {
//            double v = objects[k]->intersect(&ref_ray,rlf_color,0);
//
//            if(0<v && v<t){
//                t = v;
//                nearest = k ;
//            }
//        }
//        if(nearest!=-1){
//            t_min = objects[nearest]->intersect(&ref_ray,rlf_color,level+1);
//            for (int j=0; j<3; j++) {
//                newColor[j]+=rlf_color[j]*coEfficients[3];
//                if(newColor[j]<0) newColor[j]=0.0;
//                else if(newColor[j]>1) newColor[j]=1.0;
//            }
//        }

        return t_min;
    }
};

class Triangle : public Object{
public:
    explicit Triangle(vector<vector3D> points)
    {
        //eta = 0.1;
        for (int i = 0; i < 3; ++i) {
            tri_points[i]=points[i];
        }
    }
    void draw() override
    {
        glColor3f(color[0], color[1], color[2]);
        glBegin(GL_TRIANGLES);{
            glVertex3f(tri_points[0].x, tri_points[0].y, tri_points[0].z);
            glVertex3f( tri_points[1].x,tri_points[1].y,tri_points[1].z);
            glVertex3f(tri_points[2].x,tri_points[2].y,tri_points[2].z);

        }glEnd();
    }
    double intersect(Ray *r, double *newColor, int level) override {
        vector3D a{},b{},c{}, Ro{}, Rd{}, ab{}, ac{}, ar{};
        Ro = r->start;
        Rd = r->dir;
        a = tri_points[0];
        b = tri_points[1];
        c = tri_points[2];
        ab.set_point(a.x-b.x,a.y-b.y,a.z-b.z);
        ac.set_point(a.x-c.x,a.y-c.y,a.z-c.z);
        ar.set_point(a.x-Ro.x,a.y-Ro.y,a.z-Ro.z);

        double beta,gama,t,d,d_beta, d_gama, d_t;
        double D[3][3] = {{ab.x, ac.x, Rd.x},
                          {ab.y, ac.y, Rd.y},
                          {ab.z, ac.z, Rd.z}};

        d =      D[0][0]*(D[1][1]*D[2][2]-D[1][2]*D[2][1])
                 -D[0][1]*(D[1][0]*D[2][2]-D[2][0]*D[1][2])
                 +D[0][2]*(D[1][0]*D[2][1]-D[2][0]*D[1][1]);

        double D_beta[3][3] = {{ar.x, ac.x, Rd.x},
                               {ar.y, ac.y, Rd.y},
                               {ar.z, ac.z, Rd.z}};

        double D_gama[3][3] = {{ab.x, ar.x, Rd.x},
                               {ab.y, ar.y, Rd.y},
                               {ab.z, ar.z, Rd.z}};

        double D_t[3][3] = {{ab.x, ac.x, ar.x},
                            {ab.y, ac.y, ar.y},
                            {ab.z, ac.z, ar.z}};

        d_beta =  D_beta[0][0]*(D_beta[1][1]*D_beta[2][2]-D_beta[1][2]*D_beta[2][1])
                  -D_beta[0][1]*(D_beta[1][0]*D_beta[2][2]-D_beta[2][0]*D_beta[1][2])
                  +D_beta[0][2]*(D_beta[1][0]*D_beta[2][1]-D_beta[2][0]*D_beta[1][1]);

        d_gama =  D_gama[0][0]*(D_gama[1][1]*D_gama[2][2]-D_gama[1][2]*D_gama[2][1])
                  -D_gama[0][1]*(D_gama[1][0]*D_gama[2][2]-D_gama[2][0]*D_gama[1][2])
                  +D_gama[0][2]*(D_gama[1][0]*D_gama[2][1]-D_gama[2][0]*D_gama[1][1]);

        d_t =     D_t[0][0]*(D_t[1][1]*D_t[2][2]-D_t[1][2]*D_t[2][1])
                  -D_t[0][1]*(D_t[1][0]*D_t[2][2]-D_t[2][0]*D_t[1][2])
                  +D_t[0][2]*(D_t[1][0]*D_t[2][1]-D_t[2][0]*D_t[1][1]);

        if(abs(d)<eps)
        {
            return -1.0;
        }
        beta = d_beta/d;
        gama = d_gama/d;
        t = d_t/d;

        if(beta+gama>1 || beta<0 || gama<0 || t<0)
        {
            return -1.0;
        }
        if(level==0)
        {
            return t;
        }

        vector3D intersection_point{}, normal{};
        intersection_point = r->start + r->dir * t;
        normal = cross_product(ab,ac);
        normal.normalize();
        lighting(newColor, color, intersection_point, normal, r, level) ;

        return t;

    }
};

class Floor : public Object{
public:
    double floor_size;
    Floor(double floorWidth, double tileWidth){

        //eta = 0.3;
        ref_point.set_point(-floorWidth/2,-floorWidth/2,0.0);
        length=tileWidth;
        floor_size =floorWidth/2;
    }
    void draw() override
    {
        //cout<<floor_size<<endl;
        for(int i=(int)ref_point.x, color_switch=0; i<floor_size; i+=(int)length, color_switch=1-color_switch)
        {
            //cout<<i<<endl;
            for(int j=(int)ref_point.y, color_apply=color_switch; j<floor_size; j+=(int)length, color_apply=1-color_apply)
            {
                //cout<<j<<endl;
                glColor3f((float)color_apply,(float)color_apply,(float)color_apply);
                glBegin(GL_QUADS);{
                    glVertex3f( i, j,0);
                    glVertex3f( i+length,j,0);
                    glVertex3f(i+length,j+length,0);
                    glVertex3f(i, j+length,0);
                }glEnd();
            }
        }
        //cout<<"hi\n";
//        bool f1,f2;
//        int c1=0, c2=0;
//        for(int i=0; i<10000; i+=(int)length)
//        {
//            if(c1%2)
//            {
//                f1=true;
//                f2=false;
//            }
//            else{
//                f1=false;
//                f2=true;
//            }
//            for(int j=0; j<10000; j+=(int)length)
//            {
//                if(f1 && !f2)
//                {
//                    if(c2%2)
//                    {
//                        glColor3f(1,1,1);
//                    }
//                    else{
//                        glColor3f(1,0,0);
//                    }
//                }
//                else if(!f1 && f2)
//                {
//                    if(c2%2)
//                    {
//                        glColor3f(1,0,0);
//                    }
//                    else{
//                        glColor3f(1,1,1);
//                    }
//                }
//                glBegin(GL_QUADS);{
//                    glVertex3f( ref_point.x+i, ref_point.y+j,0);
//                    glVertex3f( ref_point.x+length+i,ref_point.y+j,0);
//                    glVertex3f(ref_point.x+length+i,ref_point.y+length+j,0);
//                    glVertex3f(ref_point.x+i, ref_point.y+length+j,0);
//                }glEnd();
//                c2++;
//            }
//            c1++;
//        }


    }

    double intersect(Ray *r, double *newColor, int level) override
    {
        vector3D normal{}, intersection_point{}, Ro{}, Rd{};
        Ro = r->start;
        Rd = r->dir;
        normal.set_point(0,0,1);

        //Ro.set_point(r->start.x-ref_point.x,r->start.y-ref_point.y,r->start.z-ref_point.z);
        //Ro.set_point(r->start.x,r->start.y,r->start.z);
        //Rd.set_point(r->dir.x,r->dir.y,r->dir.z);

        double temp,t;
        temp = dot_product(normal,Rd);
        if(abs(temp)<eps){
            return -1.0;
        }
        t = -dot_product(normal,Ro)/temp;

        if(t<0.0)
        {
            return -1.0;
        }
        intersection_point = r->start + r->dir * t;
        //intersection_point.set_point(r->start.x+r->dir.x*t,r->start.y+r->dir.y*t,r->start.z+r->dir.z*t);
        //cout<<intersection_point.x<<" "<<intersection_point.y<<" "<<intersection_point.z<<endl;

        if(intersection_point.x<-floor_size || intersection_point.x>floor_size ||
           intersection_point.y<-floor_size || intersection_point.y>floor_size)
        {
            return -1.0;
        }
        if(level==0)
        {
            return t;
        }

        double intersection_color[3];
        int i = int((intersection_point.x + floor_size) / length);
        int j = int((intersection_point.y + floor_size) / length);
        if ((i+j)%2==0) {
            intersection_color[0] = intersection_color[1] = intersection_color[2] = 0;
        }
        else {
            intersection_color[0] = intersection_color[1] = intersection_color[2] = 1;
        }
        lighting(newColor, intersection_color, intersection_point, normal, r, level) ;

        //cout<<newColor[0]<<" "<<newColor[1]<<" "<<newColor[2]<<endl;
        return t;
    }
};

void check()
{
    for(int i=0;i<lights.size();i++)
    {
        cout<<"lights: i am light\n";
    }
    for(int i=0;i<objects.size();i++)
    {
        cout<<"Objects: i am obj\n";
    }

}