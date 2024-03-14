#ifndef PARTICLESYSTEM_H
#define PARTICLESYSTEM_H

#include <vector>
#include "aParticle.h"
#include <glm.hpp>

class AParticle;
class AParticleSystem 
{
public:
    AParticleSystem();
    virtual ~AParticleSystem();


    virtual void reset();
    virtual bool isAlive(); 

	AParticle* getParticle(int index) { return mParticles[index]; }
	int getParticleNum() { return mParticles.size(); }
 
	virtual void update(double deltaT);

protected:
    std::vector<AParticle*> mParticles; // vector containing pointers to particles

public:

    bool mInfinite;                // determines whether the particle system plays once, or continuously respawns particles
    
	int mMaxParticles;    // max number of particles that can be stored in mParticles
    double mLifeSpan;              // default life span of particles
	double mDeltaT;               // size of simulation timestep
	int mSpawnTime;               // number of time steps before emmitting another particle
	int mSpawnDelay;			// count of number of time steps before emitting another particle

    glm::vec3 mStartPos;                // location of particle emitter
    glm::vec3 mStartVel;                // initial velocity of particle leaving emitter
	glm::vec3 mGravity;                 // direction and mag of gravity vector in world coords

	glm::vec3 mStartColor, mEndColor;
	double mStartAlpha, mEndAlpha;
	double mStartScale, mEndScale;
};

#endif 
