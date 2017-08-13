/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <vector>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    num_particles = 100;
    
    // init randon gen
    default_random_engine gen;
    
    // create a normal (Gaussian) distribution for x, y, theta
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    // Create particle
    Particle p;
    
    for (int i = 0; i < num_particles; ++i) {
        p.id = i;
        p.x = dist_x(gen);
        p.y = dist_y(gen);
        p.theta = dist_theta(gen);
        p.weight = 1.0; 
        // NOTE: weights may need to be at the pf level, not ind particles
        
        particles.push_back(p);
    }
    
   
   // Set init flag to true     
   is_initialized = true;
   
   // Print initial states
   for(int n=0; n<num_particles; ++n) {
       cout << "INIT: " << particles[n].id << ", x = " << particles[n].x << ", y =" << particles[n].y << ", theta = " << particles[n].theta << endl;
   }
    
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    
    // Declarations
    double x0, y0, Oo;
    
    // init randon gen
    default_random_engine gen;
    
    // create a normal (Gaussian) distribution given std_pos vector
    normal_distribution<double> dist_x(0, std_pos[0]);
    normal_distribution<double> dist_y(0, std_pos[1]);
    normal_distribution<double> dist_theta(0, std_pos[2]);
    
    // Loop through particles, and propagate states
    for (int i = 0; i < num_particles; ++i) {
               
        // Put particle states in easier to read variable names
        x0 = particles[i].x;
        y0 = particles[i].y;
        Oo = particles[i].theta; 
        
        particles[i].x = x0 + velocity/yaw_rate*(sin(Oo+yaw_rate*delta_t) - sin(Oo)) + dist_x(gen);
        particles[i].y = y0 + velocity/yaw_rate*(cos(Oo) - cos(Oo+yaw_rate*delta_t)) + dist_y(gen);
        particles[i].theta = Oo + yaw_rate*delta_t + dist_theta(gen);       
       
    }
}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    
    //Declarations
    float x_p, y_p, theta;
    float x_map, y_map;
    float range_closest, x_closest, y_closest, range_k;
    int id_closest;
    double Q, MVGaussian, MVG_mult;
    
    // For each particle, find weights  
    for (int i=0; i<particles.size(); i++) {
        
        // Init MVG_mult for weights
        MVG_mult = 1;
        
        // use local variables for particle states for ease of reading
        x_p = particles[i].x;
        y_p = particles[i].y;
        theta = particles[i].theta;
        
        // convert observations to map frame
        // based on lesson 14.13
        for (int j=0; j<observations.size(); j++) {
            
            //cout << "Obs: " << j << "  x= " << observations[j].x << "  y= " << observations[j].y << endl;
            
            // Transform observation to map frame
            //x_map = x_p + observations[j].x*cos(theta) + observations[j].y*sin(theta);
            //y_map = y_p - observations[j].x*sin(theta) + observations[j].y*cos(theta);
            
            /* from class vid https://www.youtube.com/watch?list=PLAwxTw4SYaPnfR7TzRZN-uxlxGbqxhtm2&v=-3HI3Iw3Z9g */
            x_map = x_p + (observations[j].x*cos(theta) - observations[j].y*sin(theta));
            y_map = y_p + (observations[j].x*sin(theta) + observations[j].y*cos(theta));
            
            //cout << "x_map = " << x_map << ", Y_map = " << y_map << endl;
            
            // Initialize closest range parameter with big number
            range_closest = 100.0;
            
            // find the closest landmark for the observation
            for (int k=0; k<map_landmarks.landmark_list.size(); k++) {    
           
                // Calculate distance to each landmark
                range_k = dist(x_map, y_map, map_landmarks.landmark_list[k].x_f, \
                                     map_landmarks.landmark_list[k].y_f); 
                
                
                // if current closest landmark, save 
                if (range_k < range_closest) {
                    range_closest = range_k;
                    id_closest = k+1;
                   // cout << "Replaced!  ID: " << id_closest << endl;
                    x_closest = map_landmarks.landmark_list[k].x_f;
                    y_closest = map_landmarks.landmark_list[k].y_f;  
                  
                }
            }
            //cout << "range k =" << range_k << endl;
            //cout << "range closest =" << range_closest << endl;
            //cout << "Obs ID: " << id_closest << "; x= " << x_closest << "; y= " << y_closest << "; range= " << range_closest << endl;
                
            // Use closest landmark to find MVG prob for each observation
            Q = -(pow(x_map-x_closest,2.0)/(2*pow(std_landmark[0],2.0)) + \
                        pow((y_map-y_closest),2.0)/(2*(pow(std_landmark[1],2.0))));
            MVGaussian = 1 / (2*M_PI*std_landmark[0]*std_landmark[1]) * exp(Q);
            cout << "Obs ID: " << id_closest << "; Weight = " << MVGaussian << endl;
            MVG_mult *= MVGaussian; 
           
            
        }
        cout << "Particle Weight = " << MVG_mult << endl;
        
        // update particle weight
        particles[i].weight = MVG_mult; 
        //cout << particles[i].id << ", weight = " << particles[i].weight << endl;
        
    }
    
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    // Populate weights for forming distribution
    std::vector<double> weights;
    for(int j = 0; j < num_particles; j++){
        Particle p = particles[j];
        weights.push_back(p.weight);
        //cout << particles[j].id << ", weight = " << p.weight <<endl;
    }
    
    // Initialize variables and generator
    std::vector<Particle> particles_new;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d(weights.begin(), weights.end());
    
    // Resample new set of particles
    for(int n=0; n<num_particles; ++n) {
        Particle particle_res = particles[d(gen)];
        particles_new.push_back(particle_res);    
    }
    
    // replace particles with new ones
    particles = particles_new;
    //Particle &p = particles[i];
    
    // Print updated particle states 
    for(int n=0; n<num_particles; ++n) {
        cout << "ParID: " << particles[n].id << ", x = " << particles[n].x << ", y = " << particles[n].y << ", theta = " << particles[n].theta << endl;
    }
    
 
    
    
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
