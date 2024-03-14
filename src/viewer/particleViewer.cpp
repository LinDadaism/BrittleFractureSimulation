#include "particleViewer.h"

ParticleViewer::ParticleViewer(const std::string& name) :
	Viewer(name)
{
	clearColor = ImVec4(0.6f, 0.6f, 0.6f, 1.00f);
	mParticleModel = std::make_unique<ObjModel>();
	mRocketModel = std::make_unique<ObjModel>();
	mParticleModelSphere = std::make_unique<ObjModel>();
	mParticleModel->loadObj("../obj/cube.obj");
	mRocketModel->loadObj("../obj/cone.obj");
	mParticleModelSphere->loadObj("../obj/sphere.obj");
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

ParticleViewer::~ParticleViewer()
{
}

void ParticleViewer::drawScene()
{
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	float deltaT = 1.0 / ImGui::GetIO().Framerate;
	glm::mat4 projView = mCamera.getProjView();
	drawGridGround(projView);
	mParticles.update(deltaT);
	drawParticles(projView);
	glDisable(GL_DEPTH_TEST);
}

void ParticleViewer::createGUIWindow()
{
	ImGui::Begin("Brittle Fracture Simulation Editor");
	ImGui::Text("Fracture Patterns");
	const char* fractureTypes[] = { "Uniform", "Cluster", "Brick", "Slice", "Radial", "Custom" };

	// Set spline type
	ImGui::Combo("Type", &mFractureType, fractureTypes, IM_ARRAYSIZE(fractureTypes));
	ImGui::Checkbox("Debug", &mDebug);

	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
	ImGui::End();
}

void ParticleViewer::drawParticles(const glm::mat4& projView)
{
	mModelShader->use();
	mModelShader->setMat4("uProjView", projView);
	mModelShader->setVec3("uLightPos", mLightPos);
	for (int i = 0; i < mParticles.getParticleNum(); ++i)
	{
		AParticle* particle = mParticles.getParticle(i);
		if (particle->isAlive())
		{
			// ramp the color and fade out the particle as its life runs out 
			float u = particle->getTimeToLive() / particle->m_lifeSpan;
			float scale = particle->mStartScale * u + particle->mEndScale * (1 - u);
			float alpha = particle->mStartAlpha * u + particle->mEndAlpha * (1 - u);
			glm::vec3 color = particle->mStartColor * u + particle->mEndColor * (1 - u);
			glm::vec3 pos = particle->m_Pos;
			glm::mat4 model = glm::mat4(1.0f);
			model = glm::scale(model, glm::vec3(5, 5, 5));
			model = glm::translate(model, pos);
			
			mModelShader->setMat4("uModel", model);
			mModelShader->setMat3("uModelInvTr", glm::mat3(glm::transpose(glm::inverse(model))));
			mModelShader->setVec3("color", glm::vec3(1,0,0));
			//mModelShader->setFloat("uAlpha", alpha);
			mParticleModelSphere->drawObj();
		}
	}
}