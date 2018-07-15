#include "InverseKinematics.h"


float degree2rad = 3.141592658f / 180.0f;
float rad2degree = 180.0f / 3.141592658f;

//IK workflow
void InverseKinematics::IK() {

	//check if the endeffector and the target are close enough
	Vector3 endPos = Vector3(endEffector->GlobalPos.x, endEffector->GlobalPos.y, endEffector->GlobalPos.z);
	Vector3 tarPos = Vector3(target->x, target->y, target->z);

	int i = 0;
	while ((endPos - tarPos).length() > threshold && i < iterNum)
	{
		if (mode == 0)
		{
			CCDMode();
		}
		else if (mode == 1)
		{
			JacobianMode();
			forwardKinematcis->forwardKinematicsComputation(base);
		}
		endPos = Vector3(endEffector->GlobalPos.x, endEffector->GlobalPos.y, endEffector->GlobalPos.z);
		i++;
	}
}

//CCD Mode IK for right arm
void InverseKinematics::CCDMode()
{
	//add your code here
	//hint: Do forward kinematics when the endeffector's global position need to be updated
	//with newly computed quaternions. 

	//I.e., use forwardKinematcis->forwardKinematicsComputation(base) 
	//whenever you need to update endeffector's global position.

	Joint* curr = endEffector->parent;
	//Joint* end = endEffector;
	Vector3 endPos = Vector3(endEffector->GlobalPos.x, endEffector->GlobalPos.y, endEffector->GlobalPos.z);
	Vector3 tarPos = Vector3(target->x, target->y, target->z);

	if ((endPos - tarPos).length() > threshold) {
		while (curr != base->parent) {

			Vector3 pEnd = Vector3(endEffector->GlobalPos.x, endEffector->GlobalPos.y, endEffector->GlobalPos.z);
			Vector3 pCurr = Vector3(curr->GlobalPos.x, curr->GlobalPos.y, curr->GlobalPos.z);
			Vector3 pTar = Vector3(target->x, target->y, target->z);

			Vector3 currVector = (pEnd - pCurr).normalize();
			Vector3 targetVector = (pTar - pCurr).normalize();

			float cosAngle = currVector.dot(targetVector);

			if (cosAngle > 1) {
				cosAngle = 1;
			}
			if (cosAngle < -1) {
				cosAngle = -1;
			}

			Vector3 rotVec = (currVector.cross(targetVector)).normalize();

			float angle = acos(cosAngle);
			float degree = angle * rad2degree;

			Vector4 rotQuat = forwardKinematcis->buildQuaternionRotationWithRad(angle, rotVec.x, rotVec.y, rotVec.z);
			curr->Globalquat = forwardKinematcis->quaternionMultiplication(rotQuat, curr->Globalquat);

			forwardKinematcis->forwardKinematicsComputation(curr);

			curr = curr->parent;

		}
	}


}

//Entire Right Arm Jacobian IK
void InverseKinematics::JacobianMode()
{
	//add your code here
	Vector3 pEnd = Vector3(endEffector->GlobalPos.x, endEffector->GlobalPos.y, endEffector->GlobalPos.z);
	Vector3 pTar = Vector3(target->x, target->y, target->z);

	Vector3 diffVector = (pTar - pEnd);

	// computing jacobian matrix
	Eigen::MatrixXf jac(3, 12);

	Joint* loopJoint = endEffector->parent;
	int inc = 0;
	while (loopJoint != base->parent) {
		Vector3 jointVec = Vector3(loopJoint->GlobalPos.x, loopJoint->GlobalPos.y, loopJoint->GlobalPos.z);
		Vector3 r = (pEnd - jointVec).normalize();
		Vector3 jointX = Vector3(loopJoint->GlobalXAxis.x, loopJoint->GlobalXAxis.y, loopJoint->GlobalXAxis.z);
		Vector3 jointY = Vector3(loopJoint->GlobalYAxis.x, loopJoint->GlobalYAxis.y, loopJoint->GlobalYAxis.z);
		Vector3 jointZ = Vector3(loopJoint->GlobalZAxis.x, loopJoint->GlobalZAxis.y, loopJoint->GlobalZAxis.z);

		//joint(x.cross(r)) jac matrix col
		jac(0, inc) = jointX.cross(r).x;
		jac(1, inc) = jointX.cross(r).y;
		jac(2, inc) = jointX.cross(r).z;

		//joint(y.cross(r)) jac matrix col
		jac(0, inc + 1) = jointY.cross(r).x;
		jac(1, inc + 1) = jointY.cross(r).y;
		jac(2, inc + 1) = jointY.cross(r).z;
		//joint(z.cross(r)) jac matrix col
		jac(0, inc + 2) = jointZ.cross(r).x;
		jac(1, inc + 2) = jointZ.cross(r).y;
		jac(2, inc + 2) = jointZ.cross(r).z;

		inc += 3;
		loopJoint = loopJoint->parent;
	}
	//std::cout << jac << std::endl;


	/* Computing pseudoinverse */
	Eigen::MatrixXf I = Eigen::MatrixXf::Identity(12, 12);
	double lumda = 0.1;
	Eigen::MatrixXf transJ = jac.transpose();
	Eigen::MatrixXf resultJ = ((transJ*jac + lumda * lumda*I).inverse()) * transJ;

	/* trying to get result vectors for each joint  diffVector = (pTar - pEnd);*/
	Eigen::VectorXf resVec(12, 1);
	Eigen::VectorXf pTarEndDiff(3, 1);
	pTarEndDiff << diffVector.x, diffVector.y, diffVector.z;
	resVec = resultJ * pTarEndDiff; // theta vector with 12 elements

	Joint* updateJoint = endEffector->parent;
	int index = 0;
	while (updateJoint != base->parent) {

		Vector3 resTheta = Vector3(resVec(index), resVec(index + 1), resVec(index + 2));

		Vector4 tempx, tempy, tempz;
		tempx = forwardKinematcis->buildQuaternionRotation(resTheta.x, updateJoint->GlobalXAxis.x, updateJoint->GlobalXAxis.y, updateJoint->GlobalXAxis.z);
		tempy = forwardKinematcis->buildQuaternionRotation(resTheta.y, updateJoint->GlobalYAxis.x, updateJoint->GlobalYAxis.y, updateJoint->GlobalYAxis.z);
		tempz = forwardKinematcis->buildQuaternionRotation(resTheta.z, updateJoint->GlobalZAxis.x, updateJoint->GlobalZAxis.y, updateJoint->GlobalZAxis.z);

		Vector4 quatMultiply1 = forwardKinematcis->quaternionMultiplication(tempx, tempy);
		Vector4 rotQuat = forwardKinematcis->quaternionMultiplication(quatMultiply1, tempz);
		updateJoint->Globalquat = forwardKinematcis->quaternionMultiplication(rotQuat, updateJoint->Globalquat);

		index += 3;
		updateJoint = updateJoint->parent;
	}

}