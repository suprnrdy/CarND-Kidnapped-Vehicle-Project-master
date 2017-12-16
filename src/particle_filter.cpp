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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  num_particles = 100;
  
  default_random_engine eng;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  for(int i = 0; i < num_particles; i++) {
    double sample_x, sample_y, sample_theta;
    sample_x = dist_x(eng);
    sample_y = dist_y(eng);
    sample_theta = dist_theta(eng);
    
    Particle part = {
      i,            // id
      sample_x,     // x
      sample_y,     // y
      sample_theta, // theta
      1             // weight
    };
    particles.push_back(part);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  double vdt = velocity * delta_t;
  double vyaw = velocity/yaw_rate;
  double tdt = yaw_rate * delta_t;
  
  for(int i = 0; i < num_particles; i++) {
    double sample_x, sample_y, sample_theta;
    default_random_engine eng;
    
    if(yaw_rate > 0) {
      sample_x = particles[i].x + vyaw * (sin(particles[i].theta + tdt) - sin(particles[i].theta));
      sample_y = particles[i].y + vyaw * (cos(particles[i].theta) - cos(particles[i].theta + tdt));
      sample_theta = particles[i].theta + tdt;
      normal_distribution<double> dist_theta(sample_theta, std_pos[2]);
      particles[i].theta = dist_theta(eng);
    } else {
      sample_x = particles[i].x + vdt * cos(particles[i].theta);
      sample_y = particles[i].y + vdt * sin(particles[i].theta);
    }
    
    normal_distribution<double> dist_x(sample_x, std_pos[0]);
    normal_distribution<double> dist_y(sample_y, std_pos[1]);
    
    particles[i].x = dist_x(eng);
    particles[i].y = dist_y(eng);
  }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  
  //Observations = landmarks found by lidar
  //Predicted = landmarks in the particle's orientation
  for(int i = 0; i < predicted.size(); i++) {
    double min_distance = dist(predicted[i].x, predicted[i].y, observations[0].x, observations[0].y);
    for(int l = 1; l < observations.size(); l++){
      double distance = dist(predicted[i].x, predicted[i].y, observations[i].x, observations[i].y);
      if(distance < min_distance) {
        observations[i].id = predicted[i].id;
      }
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
  //
  
  // Predict measurements to map landmarks within sensor range for each particle
  for(int p = 0; p < num_particles; p++) {
    std::vector<LandmarkObs> t_obs = observations;
    // Convert landmark observations into Map coordinate system
    for(int i = 0; i < observations.size(); i++){
      double xm = particles[i].x + (cos(-M_PI/2)*observations[i].x - (sin(-M_PI/2)*observations[i].y));
      double ym = particles[i].y + (sin(-M_PI/2)*observations[i].x + (cos(-M_PI/2)*observations[i].y));
      LandmarkObs t_ob = {observations[i].id, xm, ym};
      t_obs.push_back(t_ob);
    }
//    dataAssociation(, <#std::vector<LandmarkObs> &observations#>)
    // Calculate particle's Final weight
  }
  
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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