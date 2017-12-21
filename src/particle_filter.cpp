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
  if(!is_initialized) {
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
//      cout << "Partical Initialized: " << i << ", " << x << ", " << y << ", " << theta << endl;
      particles.push_back(part);
      weights.push_back(1);
    }
    is_initialized = true;
  }
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
    
    if(fabs(yaw_rate) > 0) {
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
//    cout << "Predicted: " << sample_x << ", " << sample_y << ", " << particles[i].theta << endl;
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  
  //Observations = landmarks found by lidar
  //Predicted = landmarks from the map
  
  
  std::vector<LandmarkObs> landmarks = predicted;
  for(int i = 0; i < observations.size(); i++) {
//    cout << "=================================" << endl;
//    cout << "Observation: " << i << endl;
    double min_distance = 1000;
    int index = 0;
    for(int l = 0; l < landmarks.size(); l++){
      double distance = dist(landmarks[l].x, landmarks[l].y, observations[i].x, observations[i].y);
      if(distance < min_distance) {
        min_distance = distance;
//        cout << landmarks[l].id << " distance: " << distance <<  endl;
        observations[i].id = landmarks[l].id;
        index = l;
      }
    }
    landmarks.erase(landmarks.begin() + index);
//    cout << "Landmark: " << observations[i].id << endl;
//    cout << "=================================" << endl;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
    std::vector<LandmarkObs> t_obs;
    std::vector<LandmarkObs> landmark_inrange;
    
    //Create list of landmarks within Sensor Range
    for(int i = 0; i < map_landmarks.landmark_list.size(); i++) {
      /*if(dist(map_landmarks.landmark_list.at(i).x_f, map_landmarks.landmark_list.at(i).y_f, particles[i].x, particles[i].y) <= sensor_range) {
        int idm = map_landmarks.landmark_list.at(i).id_i;
        double xm = map_landmarks.landmark_list.at(i).x_f;
        double ym = map_landmarks.landmark_list.at(i).y_f;
        LandmarkObs land = {idm, xm, ym};
        landmark_inrange.push_back(land);
      }*/
      int idm = map_landmarks.landmark_list.at(i).id_i;
      double xm = map_landmarks.landmark_list.at(i).x_f;
      double ym = map_landmarks.landmark_list.at(i).y_f;
      LandmarkObs land = {idm, xm, ym};
      landmark_inrange.push_back(land);
    }

//    cout << "Observations: " << observations.size() << " Landmarks: " << landmark_inrange.size() << endl;
    
    // Loop through each observation and add to t_obs if within range
    for(int i = 0; i < observations.size(); i++){
      // Convert in range observations into Map coordinate system
      double xm = particles[p].x + cos(particles[p].theta)*observations[i].x - sin(particles[p].theta)*observations[i].y;
      double ym = particles[p].y + sin(particles[p].theta)*observations[i].x + cos(particles[p].theta)*observations[i].y;
      
      LandmarkObs t_ob = {observations[i].id, xm, ym};
      t_obs.push_back(t_ob);
    }
    
    // Now do a data association to associate each observation with a landmark
    dataAssociation(landmark_inrange, t_obs);
    SetAssociations(particles[p], t_obs);
    
    // Calculate particle's Final weight
    double weight = 1;
    double rhox = std_landmark[0];
    double rhoy = std_landmark[1];
    for(int i = 0; i < t_obs.size(); i++){
      double xdiff = pow(t_obs[i].x - landmark_inrange[t_obs[i].id - 1].x, 2)/(2*pow(rhox, 2));
      double ydiff = pow(t_obs[i].y - landmark_inrange[t_obs[i].id - 1].y, 2)/(2*pow(rhoy, 2));
      double p_obs = (1/(2 * M_PI * rhox * rhoy)) * pow(M_E,-(xdiff + ydiff));
//      cout << "Observation Probability: " << p_obs << endl;
      weight *= p_obs;
      
//      cout << "Obss: " << observations[i].x << ", " << observations[i].y
//      << " t_obs: " << t_obs[i].x << ", " << t_obs[i].y << " Ass lm: " << t_obs[i].id
//      << " lm coords: " << landmark_inrange[t_obs[i].id - 1].x << ", " << landmark_inrange[t_obs[i].id - 1].y << endl;
    }
    weights[p] = weight;
    particles[p].weight = weight;
//    cout << "Particle: " << p << " weight: " << weight << endl;
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  
  // Resampling wheel
  double beta = 0.0;
  double index = int(rand() % num_particles);
  double mw = *std::max_element(weights.begin(), weights.end());

  std::vector<Particle> p_resampled;
  for(int i = 0; i < num_particles; i++) {
    beta += fmod(rand(),2 * mw);
    while(weights[index] < beta) {
      beta = beta - weights[index];
      index = fmod(index + 1, num_particles);
    }
    p_resampled.push_back(particles[index]);
  }
  particles = p_resampled;
}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<LandmarkObs>& t_obs)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
  std::vector<int> associations;
  std::vector<double> sense_x;
  std::vector<double> sense_y;
  for(int i = 0; i < t_obs.size(); i++) {
    associations.push_back(t_obs[i].id);
    sense_x.push_back(t_obs[i].x);
    sense_y.push_back(t_obs[i].y);
  }
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
