//Paul Barstad
//CS 410 Assignment 1: Transforming Objects
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <math.h>
#include "glm/glm/glm.hpp"
#include "glm/glm/gtx/string_cast.hpp"

using namespace std;

class Camera {
	public:
		int eyeX; int eyeY; int eyeZ;
		int lookX; int lookY; int lookZ;
		int upX; int upY; int upZ;
		float left; float bottom; float right; float top;
		int resWidth; int resHeight;
};
class Ambient {
	public:
		float red; float green; float blue;
};
class Constants {
	public: 
		vector<float> Ka; vector<float> Kd; vector<float> Ks;
		Ambient ambient;
};
class Light {
	public:
		int X; int Y; int Z; int w;
		float red; float green; float blue;
};
class Sphere {
	public:
		float x; float y; float z; float radius;
		float Ka_red; float Ka_green; float Ka_blue;
		float Kd_red; float Kd_green; float Kd_blue;
		float Ks_red; float Ks_green; float Ks_blue;
		float Kr_red; float Kr_green; float Kr_blue;
        float Ko_red; float Ko_green; float Ko_blue;
        float refract;
};
class Vertices {
	public:
		vector<float> ppx; vector<float> ppy; vector<float> ppz;
};
class Ray {
    public:
        float x; float y; float z;
        float vx; float vy; float vz;
};

glm::vec3 refractTRay(glm::vec3 Uinv, glm::vec3 pt, glm::vec3 N, Sphere s, float eta1, float eta2) {
    float etar = eta1 / eta2;
    float a = -1 * etar;
    float wn = dot(Uinv, N);
    float radsq = (etar * etar) * ((wn * wn) -1) + 1;
    if(radsq < 0.0) {
        //return; //T = None
    }
    else {
        float b = (etar * wn) - sqrt(radsq);
        glm::vec3 T = (a * Uinv) + (b * N);
        return T;
    }
}

Ray refractExit(glm::vec3 Uinv, glm::vec3 pt, Sphere s, float etaInside) {
    glm::vec3 center(s.x, s.y, s.z); 
    glm::vec3 norm = pt - center; norm = normalize(norm);
    glm::vec3 T1 = refractTRay(Uinv, pt, norm, s, 1.0, etaInside);
    if( (T1[0] + T1[1] + T1[2]) == 0.0) {
        //return None...
    }
    else {
        glm::vec3 exit = pt + 2 * dot((center - pt), T1) * T1;
        glm::vec3 Nin = center - exit; Nin = normalize(Nin);
        glm::vec3 T2 = refractTRay(-T1, exit, Nin, s, etaInside, 1.0);
        Ray refR;
        refR.x = exit[0]; refR.y = exit[1]; refR.z = exit[2];
        refR.vx = T2[0]; refR.vy = T2[1]; refR.vz = T2[2];
        return refR;
    }
}

glm::vec3 color(vector<Light> lightsIn, Constants consts, glm::vec3 srfNorm, glm::vec3 ptOnTri, vector<Sphere> spheres, vector<vector<int>> faces, Vertices vertices) {
	glm::vec3 ka(consts.Ka.at(0), consts.Ka.at(1), consts.Ka.at(2)); glm::vec3 amb(consts.ambient.red, consts.ambient.green, consts.ambient.blue); glm::vec3 AMBIENT = ka * amb;
	glm::vec3 diffuse(0, 0, 0); glm::vec3 testFinal;
	for(int i = 0; i < lightsIn.size(); i++) { //for each light source
		bool shadow = false; Light temp = lightsIn.at(i); glm::vec3 lpt;
        if(temp.w != 0) { lpt[0] = temp.X; lpt[1] = temp.Y; lpt[2] = temp.Z; }
        else { lpt[0] = 1000000 * temp.X; lpt[1] = 1000000 * temp.Y; lpt[2] = 1000000 * temp.Z; }
		glm::vec3 L = lpt - ptOnTri; L = normalize(L);
		//if L hits anything on the way to lpt, don't illuminate from that light
		for(int a = 0; a < spheres.size(); a++) {
			Sphere temp = spheres.at(a); float r = temp.radius;
        	glm::vec3 Cv(temp.x, temp.y, temp.z); //sphere center
        	glm::vec3 Tv = Cv - ptOnTri; 
        	float v = dot(Tv, L); float csq = dot(Tv, Tv); float disc = r*r - (csq - v*v);
        	//if disc is positive then the ray intersects the sphere
        	if(disc > 0.0) {
        		if(glm::length(Tv) > 0.00001) { //if the ray is on the circle intersecting the other side won't cause the sphere to shadow
        			if(dot(Tv, srfNorm) >= 0.0) { shadow = true; } //if dot product of srfNorm and Tv is positive light will be on the correct side
        		}
            }
		}
		/* int f = 0;
        while(shadow == false && f < faces.size()) { //don't waste time with all faces unless necessary
			for(int f = 0; f < faces.size(); f++) { //for k triangles
                int a = faces.at(f).at(0) - 1; int b = faces.at(f).at(1) -1; int c = faces.at(f).at(2) - 1;
                glm::vec3 Dv = L; glm::vec3 Lv = ptOnTri;
                //create triangle from each point
				glm::vec3 Av(vertices.ppx.at(a), vertices.ppy.at(a), vertices.ppz.at(a));// (Pprime[0][a], Pprime[1][a], Pprime[2][a]);
                glm::vec3 Bv(vertices.ppx.at(b), vertices.ppy.at(b), vertices.ppz.at(b));// (Pprime[0][b], Pprime[1][b], Pprime[2][b]);
                glm::vec3 Cv(vertices.ppx.at(c), vertices.ppy.at(c), vertices.ppz.at(c));// (Pprime[0][c], Pprime[1][c], Pprime[2][c]);
				//create matrix M for Cramer's Rule lec10n02 line 54
				glm::mat3x3 M(Av-Bv, Av-Cv, Dv); glm::mat3x3 M1(Av-Lv, Av-Cv, Dv); glm::mat3x3 M2(Av-Bv, Av-Lv, Dv); glm::mat3x3 M3(Av-Bv, Av-Cv, Av-Lv);
                //M = transpose(M); M1 = transpose(M1); M2 = transpose(M2); M3 = transpose(M3);
                if(abs(determinant(M)) > 0.00001) {
					//find beta, gamma, and t
					float beta = glm::determinant(M1)/glm::determinant(M); float gamma = glm::determinant(M2)/glm::determinant(M); float t = glm::determinant(M3)/glm::determinant(M);
					if(beta >= 0.0 && gamma >= 0.0 && (beta + gamma) <= 1.0 && t > 0.0) { shadow = true; }
				}
				f++;
			}
		} */
		float angle = glm::dot(srfNorm, L); //cos of the angle
		if(angle < 0) { /*light is on wrong side of triangle... do nothing to test*/ }
		else if(shadow == false) { //light hits face need to calculate diffuse color
			glm::vec3 kd(consts.Kd.at(0), consts.Kd.at(1), consts.Kd.at(2)); glm::vec3 I(temp.red, temp.green, temp.blue); diffuse = diffuse + (kd * I * angle);
		}
	}
	testFinal = AMBIENT + diffuse;
	for(int j = 0; j < 3; j++) { //force in bounds of 0 and 1
        if(testFinal[j] < 0.0) { testFinal[j] = 0.0; }
        if(testFinal[j] > 1.0) { testFinal[j] = 1.0; }
    }
	return testFinal * 255.0f;
}

glm::vec3 colorSPHERE(vector<Light> lightsIn, Ambient ambient, glm::vec3 srfNorm, glm::vec3 ptOnSphere, vector<Sphere> spheres, Sphere s, vector<vector<int>> faces, Vertices vertices, glm::vec3 pixpt) {
	glm::vec3 ka(s.Ka_red, s.Ka_green, s.Ka_blue); glm::vec3 kd(s.Kd_red, s.Kd_green, s.Kd_blue); glm::vec3 ks(s.Ks_red, s.Ks_green, s.Ks_blue);
	glm::vec3 amb(ambient.red, ambient.green, ambient.blue); 
	glm::vec3 AMBIENT = ka * amb; glm::vec3 diffuse(0, 0, 0); glm::vec3 specular(0, 0, 0); glm::vec3 testFinal;
    for(int i = 0; i < lightsIn.size(); i++) { //for each light source
		bool shadow = false; Light temp = lightsIn.at(i); glm::vec3 lpt;
        if(temp.w != 0) { lpt[0] = temp.X; lpt[1] = temp.Y; lpt[2] = temp.Z; }
        else { lpt[0] = 1000000 * temp.X; lpt[1] = 1000000 * temp.Y; lpt[2] = 1000000 * temp.Z; }
		glm::vec3 L = lpt - ptOnSphere; L = normalize(L);
		//if L hits anything on the way to lpt, don't illuminate from that light
		for(int a = 0; a < spheres.size(); a++) {
			glm::vec3 checkSameCenter(s.x, s.y, s.z);
			Sphere temp = spheres.at(a); float r = temp.radius;
        	glm::vec3 Cv(temp.x, temp.y, temp.z); //sphere center
        	if(checkSameCenter != Cv) {
		    	//Lv = ray starting point, Dv (Uv in sage) = ray direction
		    	glm::vec3 Tv = Cv - ptOnSphere; //base of ray to center of sphere
		    	float v = dot(Tv, L); float csq = dot(Tv, Tv); float disc = r*r - (csq - v*v);
		    	//if disc is positive then the ray intersects the sphere
		    	if(disc > 0.0) {
		    		if(glm::length(Tv) > 0.00001) { shadow = true; } //if the ray is on the circle intersecting the other side won't cause the sphere to shadow
		        }
            }
		}
		/* int f = 0;
		while(shadow == false && f < faces.size()) { //don't waste time with all faces unless necessary
			for(int f = 0; f < faces.size(); f++) { //for k triangles
                int a = faces.at(f).at(0) - 1; int b = faces.at(f).at(1) -1; int c = faces.at(f).at(2) - 1;
                //Dv used to be pixel to intersection ray and Lv used to be pixpt
                glm::vec3 Dv = L; glm::vec3 Lv = ptOnSphere;
                //create triangle from each point
				glm::vec3 Av(vertices.ppx.at(a), vertices.ppy.at(a), vertices.ppz.at(a));// (Pprime[0][a], Pprime[1][a], Pprime[2][a]);
                glm::vec3 Bv(vertices.ppx.at(b), vertices.ppy.at(b), vertices.ppz.at(b));// (Pprime[0][b], Pprime[1][b], Pprime[2][b]);
                glm::vec3 Cv(vertices.ppx.at(c), vertices.ppy.at(c), vertices.ppz.at(c));// (Pprime[0][c], Pprime[1][c], Pprime[2][c]);
				//create matrix M for Cramer's Rule lec10n02 line 54
				glm::mat3x3 M(Av-Bv, Av-Cv, Dv); glm::mat3x3 M1(Av-Lv, Av-Cv, Dv); glm::mat3x3 M2(Av-Bv, Av-Lv, Dv); glm::mat3x3 M3(Av-Bv, Av-Cv, Av-Lv);
                //M = transpose(M); M1 = transpose(M1); M2 = transpose(M2); M3 = transpose(M3);Edge
                if(abs(determinant(M)) > 0.00001) {
					//find beta, gamma, and t
					float beta = glm::determinant(M1)/glm::determinant(M);  float gamma = glm::determinant(M2)/glm::determinant(M); float t = glm::determinant(M3)/glm::determinant(M);
					if(beta >= 0.0 && gamma >= 0.0 && (beta + gamma) <= 1.0 && t > 0.0) { shadow = true; }
				}
				f++;
			}
		}*/
		float angle = glm::dot(srfNorm, L); //cos of the angle
		if(shadow == false && angle > 0.0) { //light hits face need to calculate diffuse color
			glm::vec3 I(temp.red, temp.green, temp.blue);
			diffuse = diffuse + (kd * I * angle);
            //calculate CdR
            glm::vec3 toC = pixpt - ptOnSphere; toC = normalize(toC);
            glm::vec3 spR = (2 * dot(srfNorm, L) * srfNorm) - L;
            float CdR = dot(toC, spR); float phong = pow(CdR, 20);
            if(CdR > 0.0) {
                I[0] = I[0] * phong; I[1] = I[1] * phong; I[2] = I[2] * phong;
                specular = specular + (ks * I);
            }
		}
	}
	testFinal = AMBIENT + diffuse + specular;
    testFinal[0] = testFinal[0] * s.Ko_red; testFinal[1] = testFinal[1] * s.Ko_green;
    testFinal[2] = testFinal[2] * s.Ko_blue;
	for(int j = 0; j < 3; j++) { //force in bounds of 0 and 1
        if(testFinal[j] < 0.0) { testFinal[j] = 0.0; }
        if(testFinal[j] > 1.0) { testFinal[j] = 1.0; }
    }
	return testFinal * 255.0f;
}

bool recursiveRayCast(glm::vec3 pixpt, glm::vec3 shoot, Vertices vertices, vector<vector<int>> faces, vector<Sphere> spheres, vector<Light> lights, Constants constant, int level, glm::vec3 currentColor, int v[3], string smoothVSsharp) {
	glm::vec3 Dv = shoot; //Dv represents the ray
    glm::vec3 Lv = pixpt; //Lv represents the starting point on the near plane of the ray
    float beta; float gamma; float t;
    bool hit = false; //make pixel ray hit false each time, if it hits it will continue to search for the lowest t (closest Face) 
    bool closerSphere = false; //also check for closer/general hits on spheres
    float smallestT = 1000000.0; int closestFace; glm::vec3 surfNorm; glm::vec3 ptInTri; 
    int bestA; int bestB; int bestC; float bestBeta; float bestGamma;
    for(int f = 0; f < faces.size(); f++) { //for k triangles
        int a = faces.at(f).at(0) - 1; int b = faces.at(f).at(1) -1; int c = faces.at(f).at(2) - 1;
        //create triangle from each point
		glm::vec3 Av(vertices.ppx.at(a), vertices.ppy.at(a), vertices.ppz.at(a));// (Pprime[0][a], Pprime[1][a], Pprime[2][a]);
        glm::vec3 Bv(vertices.ppx.at(b), vertices.ppy.at(b), vertices.ppz.at(b));// (Pprime[0][b], Pprime[1][b], Pprime[2][b]);
		glm::vec3 Cv(vertices.ppx.at(c), vertices.ppy.at(c), vertices.ppz.at(c));// (Pprime[0][c], Pprime[1][c], Pprime[2][c]);
		//create matrix M for Cramer's Rule lec10n02 line 54
		glm::mat3x3 M(Av-Bv, Av-Cv, Dv); glm::mat3x3 M1(Av-Lv, Av-Cv, Dv); glm::mat3x3 M2(Av-Bv, Av-Lv, Dv); glm::mat3x3 M3(Av-Bv, Av-Cv, Av-Lv);
        //M = transpose(M); M1 = transpose(M1); M2 = transpose(M2); M3 = transpose(M3);
        if(abs(determinant(M)) > 0.00001) {
			//find beta, gamma, and t
			beta = glm::determinant(M1)/glm::determinant(M); gamma = glm::determinant(M2)/glm::determinant(M); t = glm::determinant(M3)/glm::determinant(M);
			if(beta >= 0.0 && gamma >= 0.0 && (beta + gamma) <= 1.0 && t > 0.0) {
                //cout << "recursive face hit" << endl;
				glm::vec3 AB = Bv - Av; glm::vec3 AC = Cv - Av; glm::vec3 norm = cross(AB, AC); norm = normalize(norm);	
				float dt = dot(Dv, norm);
				if(dt > 0.0) { norm = -norm; }	//if dot product less than zero flip it		
				if(hit == false) { //new face hit 
					hit = true; closestFace = f; smallestT = t; surfNorm = norm; ptInTri = Lv + Dv*t;
                    bestA = a; bestB = b; bestC = c; bestBeta = beta; bestGamma = gamma;
				}
				else if(smallestT > t) { //closer face
					smallestT = t; closestFace = f; surfNorm = norm; ptInTri = Lv + Dv*t;
                    bestA = a; bestB = b; bestC = c; bestBeta = beta; bestGamma = gamma;
				}
			}
		} //end first determinant of M check
    }
    Sphere closest;
    glm::vec3 closestSN; glm::vec3 closestPT;
    for(int s = 0; s < spheres.size(); s++) { //for each sphere
    	//check if sphere intersection
    	Sphere temp = spheres.at(s); float r = temp.radius;
    	glm::vec3 Cv(temp.x, temp.y, temp.z); //sphere center
    	//Lv = ray starting point, Dv (Uv in sage) = ray direction
    	glm::vec3 Tv = Cv - Lv; //base of ray to center of sphere
    	float v = dot(Tv, Dv); float csq = dot(Tv, Tv); float disc = r*r - (csq - v*v); float d = sqrt(disc);
    	if(disc > 0.0 && (v - d) > 0.000001) {//if disc is positive then the ray intersects the sphere
            //check for same sphere intersection with v-d above
    		glm::vec3 pt = Lv + (v-d) * Dv; glm::vec3 sn = pt - Cv; sn = normalize(sn);
    		glm::vec3 toSphere = pt - Lv; float testT = glm::length(toSphere);
    		//need to check if closer than other faces, if not color sphere...
	    	if(hit == true ) {
	    		if(testT < smallestT) { //sphere is closer than face
	    			closerSphere = true; smallestT = testT; closest = temp; closestSN = sn; closestPT = pt;
				}
			}
			else if(testT < smallestT) { //no face intersection
                closerSphere = true; smallestT = testT; closest = temp; closestSN = sn; closestPT = pt;
			}
        }
    }
    if(hit == true && closerSphere == false && smoothVSsharp == "smooth") {
        //cout << "smoothing " << surfNorm[0] << " " << surfNorm[1] << " " << surfNorm[2] << endl;
        //smooth normals
        vector<vector<int>> shareA; vector<vector<int>> shareB; vector<vector<int>> shareC;
        for(int i = 0; i < faces.size(); i++) {
            //if face vertex a matches add it to shareA vector
            if(faces.at(i).at(0) == bestA || faces.at(i).at(1) == bestA || faces.at(i).at(2) == bestA) { shareA.push_back(faces.at(i)); }
            //if face vertex b matches add it to shareB vector
            if(faces.at(i).at(0) == bestB || faces.at(i).at(1) == bestB || faces.at(i).at(2) == bestB) { shareB.push_back(faces.at(i)); }
            //if face vertex c matches add it to shareC vector
            if(faces.at(i).at(0) == bestC || faces.at(i).at(1) == bestC || faces.at(i).at(2) == bestC) { shareC.push_back(faces.at(i)); }
        }
        int countA = 0;
        //for each face in shareA
        glm::vec3 temp = surfNorm;
        for(int i = 0; i < shareA.size(); i++) {
            //calculate Ni = cross(A-B, A-C) and normalize it
            int a = shareA.at(i).at(0) - 1; int b = shareA.at(i).at(1) -1; int c = shareA.at(i).at(2) - 1;
            glm::vec3 Av(vertices.ppx.at(a), vertices.ppy.at(a), vertices.ppz.at(a));
            glm::vec3 Bv(vertices.ppx.at(b), vertices.ppy.at(b), vertices.ppz.at(b));
            glm::vec3 Cv(vertices.ppx.at(c), vertices.ppy.at(c), vertices.ppz.at(c));
            glm::vec3 AB = Bv - Av; glm::vec3 AC = Cv - Av; glm::vec3 Ni = cross(AB, AC); Ni = normalize(Ni);
            //add Ni to surfNorm
            if(dot(surfNorm, Ni) > 0.3927) {
                //cout << "found a" << endl;
                temp = temp + Ni;
                //increment countA
                countA++;
            }
        }
        glm::vec3 Na;
        if(countA > 0) { Na = temp; Na[0] = Na[0]/countA; Na[1] = Na[1]/countA; Na[2] = Na[2]/countA; }
        else { Na = surfNorm; }
        
        int countB = 0; temp = surfNorm;
        for(int i = 0; i < shareB.size(); i++) {
            //calculate Ni = cross(A-B, A-C) and normalize it
            int a = shareB.at(i).at(0) - 1; int b = shareB.at(i).at(1) -1; int c = shareB.at(i).at(2) - 1;
            glm::vec3 Av(vertices.ppx.at(a), vertices.ppy.at(a), vertices.ppz.at(a));
            glm::vec3 Bv(vertices.ppx.at(b), vertices.ppy.at(b), vertices.ppz.at(b));
            glm::vec3 Cv(vertices.ppx.at(c), vertices.ppy.at(c), vertices.ppz.at(c));
            glm::vec3 AB = Bv - Av; glm::vec3 AC = Cv - Av; glm::vec3 Ni = cross(AB, AC); Ni = normalize(Ni);
            //add Ni to surfNorm
            if(dot(surfNorm, Ni) > 0.3927) {
                //cout << "found b" << endl;
                temp = temp + Ni;
                //increment countA
                countB++;
            }
        }
        glm::vec3 Nb;
        if(countB > 0) { Nb = temp; Nb[0] = Nb[0]/countB; Nb[1] = Nb[1]/countB; Nb[2] = Nb[2]/countB; }
        else { Nb = surfNorm; }
        
        int countC = 0; temp = surfNorm;
        for(int i = 0; i < shareC.size(); i++) {
            //calculate Ni = cross(A-B, A-C) and normalize it
            int a = shareC.at(i).at(0) - 1; int b = shareC.at(i).at(1) -1; int c = shareC.at(i).at(2) - 1;
            glm::vec3 Av(vertices.ppx.at(a), vertices.ppy.at(a), vertices.ppz.at(a));
            glm::vec3 Bv(vertices.ppx.at(b), vertices.ppy.at(b), vertices.ppz.at(b));
            glm::vec3 Cv(vertices.ppx.at(c), vertices.ppy.at(c), vertices.ppz.at(c));
            glm::vec3 AB = Bv - Av; glm::vec3 AC = Cv - Av; glm::vec3 Ni = cross(AB, AC); Ni = normalize(Ni);
            //add Ni to surfNorm
            if(dot(surfNorm, Ni) > 0.3927) {
                //cout << "found c" << endl;
                temp = temp + Ni;
                //increment countA
                countC++;
            }
        }
        glm::vec3 Nc;
        if(countC > 0) { Nc = temp; Nc[0] = Nc[0]/countC; Nc[1] = Nc[1]/countC; Nc[2] = Nc[2]/countC; }
        else { Nc = surfNorm; }
        
        //set surfNorm to normInterpolated
        float aChange = 1 - bestBeta - bestGamma; Na[0] = Na[0]*aChange; Na[1] = Na[1]*aChange; Na[2] = Na[2]*aChange;
        Nb[0] = Nb[0]*bestBeta; Nb[1] = Nb[1]*bestBeta; Nb[2] = Nb[2]*bestBeta;
        Nc[0] = Nc[0]*bestGamma; Nc[1] = Nc[1]*bestGamma; Nc[2] = Nc[2]*bestGamma;
        
        surfNorm = Na + Nb + Nc;
        //cout << "new " << surfNorm[0] << " " << surfNorm[1] << " " << surfNorm[2] << endl;
    }
    
    glm::vec3 newpt;
    glm::vec3 refR;
    if(closerSphere == true) { //color sphere
        Sphere temp = closest;                        
        currentColor += colorSPHERE(lights, constant.ambient, closestSN, closestPT, spheres, temp, faces, vertices, pixpt);
        newpt = closestPT;
        glm::vec3 Uinv = -shoot;
        refR = 2 * dot(closestSN, Uinv) * closestSN - Uinv; refR = normalize(refR);
    }
    if(hit == true && closerSphere == false) { //color face
        //cout << "color face" << endl;
    	currentColor += color(lights, constant, surfNorm, ptInTri, spheres, faces, vertices);
    	newpt = ptInTri; 
        glm::vec3 Uinv = -shoot;
    	refR = 2 * dot(surfNorm, Uinv) * surfNorm - Uinv; refR = normalize(refR);
    }
	else if(closerSphere == false) { //no hit with ray
        v[0] = currentColor[0]; v[1] = currentColor[1]; v[2] = currentColor[2];
        return false;
	}
	
	
	if(level > 0 && currentColor[0] + currentColor[1] + currentColor[2] > 0) { //need to recurse again
        recursiveRayCast(newpt, refR, vertices, faces, spheres, lights, constant, level-1, currentColor, v, smoothVSsharp);
    }
    //need if check for closerSphere in case using models??
    if(level > 0 && closest.Ko_red + closest.Ko_green + closest.Ko_blue < 3.0 && currentColor[0] + currentColor[1] + currentColor[2] > 0) {
        //glm::vec3 thru(0, 0, 0);
        Ray r = refractExit(-shoot, newpt, closest, closest.refract);
        glm::vec3 fraRpt(r.x, r.y, r.z);
        glm::vec3 fraR(r.vx, r.vy, r.vz);
        recursiveRayCast(fraRpt, fraR, vertices, faces, spheres, lights, constant, level-1, currentColor, v, smoothVSsharp);
        //currentColor[0] = currentColor[0] * (1.0 - closest.Ko_red);
        //currentColor[1] = currentColor[1] * (1.0 - closest.Ko_green);
        //currentColor[2] = currentColor[2] * (1.0 - closest.Ko_blue);
    }
	if(level == 0) { //done recursing
		int red = currentColor[0]; int green = currentColor[1]; int blue = currentColor[2];
		if(red > 255) { red = 255; } if(green > 255) { green = 255; } if(blue > 255) { blue = 255; }
		v[0] = red; v[1] = green; v[2] = blue;
        return true;
	}
}


int main(int argc, char *argv[]) {
	//read driver file in
	fstream driver; driver.open(argv[1]);
    string outName = argv[2]; string line;
    int dVal; int recursionLevel;  
    float wx; float wy; float wz;
    float theta; float scale;
    float tx; float ty; float tz;
    string smoothVSsharp; string modelFile;
    Camera camera = Camera(); Ambient ambient = Ambient(); vector<Light> lights; vector<Sphere> spheres;
    while(getline(driver, line)) { 
        //save all variables into strings
        if(line.length() != 0) {
		    if(line.substr(0, 3) == "eye") { //if eye
		        vector<int> eye; int lastSpace = 4;
		        for(unsigned int i = 4; i < line.length(); i++) { //start at 4 to skip "eye "
		            if(line.at(i) == ' ') { eye.push_back(stoi(line.substr(lastSpace, i-lastSpace))); lastSpace = i+1; }
		            if(i == line.length()-1) { eye.push_back(stoi(line.substr(lastSpace, i-lastSpace+1))); }
		        }
		        camera.eyeX = eye.at(0); camera.eyeY = eye.at(1); camera.eyeZ = eye.at(2);
		    }
		    if(line.substr(0, 4) == "look") { //if look
		        vector<int> look; int lastSpace = 5;
		        for(unsigned int i = 5; i < line.length(); i++) { //start at 5 to skip "look "
		            if(line.at(i) == ' ') { look.push_back(stoi(line.substr(lastSpace, i-lastSpace))); lastSpace = i+1; }
		            if(i == line.length()-1) { look.push_back(stoi(line.substr(lastSpace, i-lastSpace+1))); }
		        }
		        camera.lookX = look.at(0); camera.lookY = look.at(1); camera.lookZ = look.at(2);
		    }
		    if(line.substr(0, 2) == "up") {//if up
		        vector<int> up; int lastSpace = 3;
				for(unsigned int i = 3; i < line.length(); i++) { //start at 3 to skip "up "
		            if(line.at(i) == ' ') { up.push_back(stoi(line.substr(lastSpace, i-lastSpace))); lastSpace = i+1; }
		            if(i == line.length()-1) { up.push_back(stoi(line.substr(lastSpace, i-lastSpace+1))); }
		        }
		        camera.upX = up.at(0); camera.upY = up.at(1); camera.upZ = up.at(2);
		    }
		    if(line.substr(0, 1) == "d") {//if d
		        dVal = stoi(line.substr(2, line.length()-2));
		    }
		    if(line.substr(0, 6) == "bounds") {//if bounds
		        vector<float> bounds; int lastSpace = 7;
				for(unsigned int i = 7; i < line.length(); i++) { //start at 7 to skip "bounds "
		            if(line.at(i) == ' ') { bounds.push_back(stof(line.substr(lastSpace, i-lastSpace))); lastSpace = i+1; }
		            if(i == line.length()-1) { bounds.push_back(stof(line.substr(lastSpace, i-lastSpace+1))); }
		        }
                camera.left = bounds.at(0); camera.bottom = bounds.at(1); camera.right = bounds.at(2); camera.top = bounds.at(3);
		    }
		    if(line.substr(0, 3) == "res") {//if res
		        vector<int> res; int lastSpace = 4;
		        for(unsigned int i = 4; i < line.length(); i++) { //start at 4 to skip "res "
		            if(line.at(i) == ' ') { res.push_back(stoi(line.substr(lastSpace, i-lastSpace))); lastSpace = i+1; }
		            if(i == line.length()-1) { res.push_back(stoi(line.substr(lastSpace, i-lastSpace+1))); }
		        }
		        camera.resWidth = res.at(0); camera.resHeight = res.at(1);
		    }
		    if(line.substr(0, 14) == "recursionLevel") {
                recursionLevel = stoi(line.substr(15, line.length()-15));
            }
		    if(line.substr(0, 7) == "ambient") {//if ambient
		        vector<float> amb; int lastSpace = 8;
		        for(unsigned int i = 8; i < line.length(); i++) { //start at 8 to skip "ambient "
		            if(line.at(i) == ' ') { amb.push_back(stod(line.substr(lastSpace, i-lastSpace))); lastSpace = i+1; }
		            if(i == line.length()-1) { amb.push_back(stod(line.substr(lastSpace, i-lastSpace+1))); }
		        } 
		        ambient.red = amb.at(0); ambient.green = amb.at(1); ambient.blue = amb.at(2);
		    }
		    if(line.substr(0, 5) == "light") {//if light (could be multiple)
		        //if light.size() > 0 then multiple lights in file...
		        vector<float> light; int lastSpace = 6;
		        for(unsigned int i = 6; i < line.length(); i++) { //start at 6 to skip "light "
		            if(line.at(i) == ' ') { light.push_back(stod(line.substr(lastSpace, i-lastSpace))); lastSpace = i+1; }
		            if(i == line.length()-1) { light.push_back(stod(line.substr(lastSpace, i-lastSpace+1))); }
		        }
		        Light l = Light();
		        for(unsigned int i = 0; i < light.size(); i++) { //check if problem for more than one light
		            if(i%7 == 0) { l.X = light.at(i); } if(i%7 == 1) { l.Y = light.at(i); } if(i%7 == 2) { l.Z = light.at(i); }
		        	if(i%7 == 3) { l.w = light.at(i); } 
		        	if(i%7 == 4) { l.red = light.at(i); } if(i%7 == 5) { l.green = light.at(i); } if(i%7 == 6) { l.blue = light.at(i); lights.push_back(l); }
		        }
		    }
		    if(line.substr(0, 5) == "model") { //if model
		        vector<string> results; int lastSpace = 6;
		        for(unsigned int i = 6; i < line.length(); i++) { //start at 6 to skip "model "
		            if(line.at(i) == ' ') { results.push_back(line.substr(lastSpace, i-lastSpace)); lastSpace = i+1; }
		        }
		        results.push_back(line.substr(lastSpace, line.length()-1)); // get the obj file name
		        wx = stod(results.at(0)); wy = stod(results.at(1)); wz = stod(results.at(2));
		        theta = stod(results.at(3)); scale = stod(results.at(4));
		        tx = stod(results.at(5)); ty = stod(results.at(6)); tz = stod(results.at(7));
                smoothVSsharp = results.at(8); modelFile = results.at(9);
		    }
			if(line.substr(0, 6) == "sphere") {//if sphere
		        vector<float> sphereInfo; int lastSpace = 7;
				for(unsigned int i = 7; i < line.length(); i++) { //start at 7 to skip "sphere "
		            if(line.at(i) == ' ') { sphereInfo.push_back(stof(line.substr(lastSpace, i-lastSpace))); lastSpace = i+1; }
		            if(i == line.length()-1) { sphereInfo.push_back(stof(line.substr(lastSpace, i-lastSpace+1))); }
		        }
				Sphere tempSphere = Sphere();
				tempSphere.x = sphereInfo.at(0); tempSphere.y = sphereInfo.at(1); tempSphere.z = sphereInfo.at(2); 
				tempSphere.radius = sphereInfo.at(3);
				tempSphere.Ka_red = sphereInfo.at(4); tempSphere.Ka_green = sphereInfo.at(5); tempSphere.Ka_blue = sphereInfo.at(6); 
				tempSphere.Kd_red = sphereInfo.at(7); tempSphere.Kd_green = sphereInfo.at(8); tempSphere.Kd_blue = sphereInfo.at(9);
				tempSphere.Ks_red = sphereInfo.at(10); tempSphere.Ks_green = sphereInfo.at(11); tempSphere.Ks_blue = sphereInfo.at(12); 
				tempSphere.Kr_red = sphereInfo.at(13); tempSphere.Ks_green = sphereInfo.at(14); tempSphere.Ks_blue = sphereInfo.at(15);
                tempSphere.Ko_red = sphereInfo.at(16); tempSphere.Ko_green = sphereInfo.at(17); tempSphere.Ko_blue = sphereInfo.at(18);
                tempSphere.refract = sphereInfo.at(19);
				spheres.push_back(tempSphere); //sphereToString(tempSphere);
			}
		}
    }
    //------------------------------------------INPUT COMPLETE--------------------------------------
    
    
        //create MM matrix
        glm::vec3 W(wx, wy, wz); W = normalize(W); //make W unit length
        //create M vector that is not parallel to W
        glm::vec3 M(W[0], W[1], W[2]); //constant
        float Mmin; int spot = 0;
        for(int z = 0; z < 3; z++) {
            if(z == 0) { Mmin = M[0]; }
            if(M[z] <= Mmin) { Mmin = M[z]; spot = z; }
        }
        M[spot] = 1.0;
        glm::vec3 U = cross(W, M); U = normalize(U); //constant
        glm::vec3 V = cross(W, U); //constant
        
        glm::mat3 RM(U[0], U[1], U[2], V[0], V[1], V[2], W[0], W[1], W[2]); 
        glm::mat3 RMt = transpose(RM); //RMt is constant
        glm::mat3 RMRMt = RM * RMt; //identity matrix
        
        float arad = (theta / 180) * M_PI; float ca = cos(arad); float sa = sin(arad);
        glm::mat3 RMz = RMRMt; 
        RMz[0][0] = ca; RMz[0][1] = -1*sa; RMz[1][0] = sa; RMz[1][1] = ca;
        
        glm::mat3 RT = RM * RMz * RMt;//RMt * RMz * RM; commented out and reversed the order which produced the same matrix as jupyter notebook
        
        //Matrix multiply against .obj file...
        fstream obj; obj.open(modelFile);
        string current; string topComments; string mtlFile;
        vector<string> verticies; vector<string> fLines;
        //sort original .obj file
        while(getline(obj, current)) { 
        	if(current.length() != 0) {
		        if(current[0] == '#') { topComments += current + '\n'; } //if comments at top of file store them
		        else if(current.substr(0,6) == "mtllib") { mtlFile = current.substr(7, current.length()-7); }//store the material details file name 
		        else if(current[0] == 'v' && current[1] != 'n' && current[1] != 't') { verticies.push_back(current.substr(2, current.length()-1)); } //if line starts with v (not vn) store the vertex
		        else if(current[0] == 'f') { fLines.push_back(current + '\n'); }
			}
        }
        
        fstream mtl; mtl.open(mtlFile);
        string ka; string kd; string ks;
        while(getline(mtl, current)) {
        	if(current.length() != 0) {
		    	if(current.substr(0,2) == "Ka") { ka = current; }
		    	if(current.substr(0,2) == "Kd") { kd = current; }
		    	if(current.substr(0,2) == "Ks") { ks = current; }
		    }
        }
        
        vector<float> Ka; //filling with ambient values of material
        int lastSpace = 3;
        for(unsigned int i = 3; i < ka.length(); i++) { //start at 3 to skip "Ka "
        	if(ka.at(i) == ' ') {
        	    Ka.push_back(stod(ka.substr(lastSpace, i-lastSpace)));
                lastSpace = i+1;
            }
            if(i == ka.length()-1) {
                Ka.push_back(stod(ka.substr(lastSpace, i-lastSpace+1)));
            }
        }
        vector<float> Kd; //filling with diffuse values of material
        lastSpace = 3;
        for(unsigned int i = 3; i < kd.length(); i++) { //start at 3 to skip "Kd "
        	if(kd.at(i) == ' ') {
        	    Kd.push_back(stod(kd.substr(lastSpace, i-lastSpace)));
                lastSpace = i+1;
            }
            if(i == kd.length()-1) {
                Kd.push_back(stod(kd.substr(lastSpace, i-lastSpace+1)));
            }
        }
       vector<float> Ks; //filling with specular values of material
       lastSpace = 3;
       for(unsigned int i = 3; i < ks.length(); i++) {
       		if(ks.at(i) == ' ') { Ks.push_back(stod(ks.substr(lastSpace, i-lastSpace))); lastSpace = i+1; }
            if(i == ks.length()-1) { Ks.push_back(stod(ks.substr(lastSpace, i-lastSpace+1))); }
       }
       Constants constants;
       constants.Ka = Ka; constants.Kd = Kd; constants.Ks = Ks; constants.ambient = ambient;
        
        vector<vector<int>> faces;
        for(unsigned int i = 0; i < fLines.size(); i++) { //loop through f lines which contain triangles of points
        	//cout << fLines.at(i) << endl;
            int j = 2;
            string temp1 = "";
            if(j == 2) { while(fLines.at(i).at(j) != '/') { j++; } temp1 = fLines.at(i).substr(2,j-2); }
            //find whitespace + 1 until '//'
            while(fLines.at(i).at(j) != ' ') { j++;}
            j++; int space = j; //skip whitespace 
            while(fLines.at(i).at(j) != '/') { j++; } string temp2 = fLines.at(i).substr(space,j-space);
            //find another whitespace + 1 until '//'
            while(fLines.at(i).at(j) != ' ') { j++; }
            j++; space = j; //skip whitespace
            while(fLines.at(i).at(j) != '/') { j++; } string temp3 = fLines.at(i).substr(space,j-space);
            vector<int> temp;
        	temp.push_back(stoi(temp1)); temp.push_back(stoi(temp2)); temp.push_back(stoi(temp3));
        	faces.push_back(temp); //each index in faces has a three vertexes, making each triangle to loop through
        }
        
        vector<float> points;
        for(unsigned int d = 0; d < verticies.size(); d++) { //loop through strings of verticies (each string is one vertex)
            int lastSpace = 0;
            string temp = verticies.at(d);
            //cout << temp << endl;
            for(int e = 0; e < temp.length(); e++) { 
                if(temp[e] == ' ') { //if you find a space
                    string point = temp.substr(lastSpace, e-lastSpace);
                    lastSpace = e+1;
                    points.push_back(stod(point)); //put the x or y value into the points vector
                }
                if(e == temp.length()-1) { //if you reach the end of the string
                    string point = temp.substr(lastSpace, e-lastSpace);
                    points.push_back(stod(point)); //put the z value into the points vector 
                }
            }
        }
        
        spot = 0;
        float P[3][points.size()/3]; //make a matrix that is 3x(number of verticies)
        int rowCount = 0;
        int columnCount = 0;
        for(unsigned int f = 0; f < points.size(); f++) {
            if(f % 3 == 0) { P[rowCount][columnCount] = points.at(f); rowCount++; } //x value
            if(f % 3 == 1) { P[rowCount][columnCount] = points.at(f); rowCount++; } //y value
            if(f % 3 == 2) { P[rowCount][columnCount] = points.at(f); rowCount = 0; columnCount++; } //z value
        } //P now contains the original verticies...
        
        float Pprime[3][points.size()/3];
        for(int a = 0; a < 3; a++) {
        	for(int b = 0; b < points.size()/3; b++) { //initialize each spot in Pprime to zero
        		Pprime[a][b] = 0.0;
        	}
        }
        for(int y = 0; y < 3; y++) { //for each row in RT
            for(int z = 0; z < points.size()/3; z++) { //for each spot in row of RT = each spot in column of points
                for(int r = 0; r < 3; r++) { //for each column in points
                    Pprime[y][z] += RT[y][r] * P[r][z]; //do the multiplication against each value
                    if(Pprime[y][z] < 0.01 && Pprime[y][z] > -0.01) {
                        Pprime[y][z] = 0.0;
                    }
                }
            }
        }
        //scale object
        for(int g = 0; g < 3; g++) {
            for(int h = 0; h < points.size()/3; h++) {
                Pprime[g][h] = Pprime[g][h] * scale;
            } 
        }
        //translate object
        for(int j = 0; j < points.size()/3; j++) {
            Pprime[0][j] = Pprime[0][j] + tx; Pprime[1][j] = Pprime[1][j] + ty; Pprime[2][j] = Pprime[2][j] + tz;
        }
        //-------------------------------------------Transformation of object complete-------------------
        
        
        vector<float> ppx; vector<float> ppy; vector<float> ppz;
        for(int d = 0; d < points.size()/3; d++) {
        	ppx.push_back(Pprime[0][d]); ppy.push_back(Pprime[1][d]); ppz.push_back(Pprime[2][d]);
        }
        Vertices vertices = Vertices(); vertices.ppx = ppx; vertices.ppy = ppy; vertices.ppz = ppz;      	
        	
        glm::vec3 EV(camera.eyeX, camera.eyeY, camera.eyeZ);	
        glm::vec3 LV(camera.lookX, camera.lookY, camera.lookZ);
        glm::vec3 WV = EV - LV; WV = normalize(WV);
        glm::vec3 UP(camera.upX, camera.upY, camera.upZ);
        glm::vec3 UV = cross(UP, WV); UV = normalize(UV);
        glm::vec3 VV = cross(WV, UV); VV = normalize(VV);
        
        ofstream outFile; outFile.open(outName);
        outFile << "P3\n" << camera.resWidth << " " << camera.resHeight << " 255" << endl; 
        
        int hits = 0; float near = -1.0*dVal; 
        for(float j = 0; j < camera.resWidth; j++) { //for pixels across
            for(float i = 0; i < camera.resHeight; i++) { //for pixels down
                float px = i/(camera.resWidth - 1) * (camera.right - camera.left) + camera.left;
                float py = j/(camera.resHeight - 1) * (camera.bottom - camera.top) + camera.top;
                //-----------------------------------------Raycast-------------------------------------------
                glm::vec3 pixpt = EV + (near * WV) + (px * UV) + (py * VV); 
                glm::vec3 shoot = pixpt - EV; shoot = normalize(shoot);
                
                int v[3];
                glm::vec3 currentColor(0, 0, 0);
                bool recursiveHit = false;
                recursiveHit = recursiveRayCast(pixpt, shoot, vertices, faces, spheres, lights, constants, recursionLevel, currentColor, v, smoothVSsharp);
                outFile << v[0] << " " << v[1] << " " << v[2] << " ";
                
            } 
            outFile << endl;   
        }
        
	return 0;
}


// sphere 35 60 20 9 1.0 0.0 0.0 1.0 0.0 0.0 1.0 1.0 1.0 0.9 0.9 0.9 1.0 1.0 1.0 1.3
// sphere 50 35 20 9 0.0 0.0 1.0 0.0 0.0 1.0 1.0 1.0 1.0 0.9 0.9 0.9 1.0 1.0 1.0 1.3
// sphere 65 60 20 9 0.0 1.0 0.0 0.0 1.0 0.0 1.0 1.0 1.0 0.9 0.9 0.9 1.0 1.0 1.0 1.3
// sphere 50 50 50 9 0.2 0.2 0.2 0.6 0.6 0.6 0.5 0.5 0.5 0.9 0.9 0.9 0.2 0.2 0.2 2.0
