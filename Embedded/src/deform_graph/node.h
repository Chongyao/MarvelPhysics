#ifndef _NODE_H_
#define _NODE_H_

#include <set>
#include <Eigen/Dense>

using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Matrix3d;

class Node
{
public:

    // Default constructor.
    Node();

    Node(Vector3d _position);

    ~Node();

    void addDeltaRotation(Matrix3d &delta);

    void addDeltaTranslation(Vector3d &delta);

    Vector3d getVelocityFrame() const;

    void setVelocityFrame(Vector3d velocity);

    // void setPosition(Vector3d _position);

    Vector3d getPosition() const;

    Matrix3d matRotation() const;

    Vector3d getTranslation() const;

    Vector3d getTranslationFrame() const;

    Vector3d updateTranslationFrame();

    Vector3d applyMapping(Vector3d &p);

    void setTransformation(Matrix3d &_rotation, Vector3d &_translation);

    Vector3d transformPosition(Vector3d &vpos);

    Vector3d transformNormal(Vector3d &normal);

    void addNeighbor(Node * n);

    std::set<Node *> getNeighbors();

    // Rotation term: not currently used
    double getRotValue();

    // Regularization term: not currently used
    double getRegValue();

    // Get [(c1*c2) (c1*c3) (c2*c3) (c1*c1-1) (c2*c2-1) (c3*c3-1)]
    VectorXd getRotTerm();

    // Get all neighbor's [Rj * (gk - gj) + gj + tj - (gk + tk)]
    MatrixXd getRegTerm();

    // Get a certain neighbor's [Rj * (gk - gj) + gj + tj - (gk + tk)]
    Vector3d getRegTerm(Node * neighbor);

private:
	bool transformed; // if the node has been transformed
	const Vector3d position;
	Matrix3d rotation;
	Vector3d translation;
    Vector3d translation_f; // translation per frame
	std::set<Node *> neighbors; // Neighbor nodes
    Vector3d velocity_f; // used for animation
};

#endif
