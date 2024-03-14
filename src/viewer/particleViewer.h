#pragma once
#include "viewer.h"
#include "aParticleSystem.h"
#include "objmodel.h"

class ParticleViewer : public Viewer
{
public:
	ParticleViewer(const std::string& name);
	virtual ~ParticleViewer();

	virtual void drawScene() override;
	virtual void createGUIWindow() override;

private:
	int mDemo = 0;	// 0 for particles, 1 for fireworks
	int mParticleModelType = 0; // 0 for cube, 1 for sphere

	int mFractureType;
	bool mDebug = true;

	AParticleSystem mParticles;

	std::unique_ptr<ObjModel> mParticleModel;
	std::unique_ptr<ObjModel> mParticleModelSphere;
	std::unique_ptr<ObjModel> mRocketModel;

	void drawParticles(const glm::mat4& projView);
};