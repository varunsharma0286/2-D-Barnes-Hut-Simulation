#include<conio.h>
#include<stdio.h>
#include <time.h>
#include <vector>
#include<iostream>
using namespace std;

#define THETA 0.5
#define XDIM 16
#define YDIM 16
#define ZDIM 1024
#define XMIN 0.0
#define YMIN 0.0
#define MASS 1.0
#define INITVELOCITY 0.0
#define START 0
#define END 1

#define NW 1
#define NE 2
#define SW 3
#define SE 4
#define INVALID -1

#define G 6.673e-11
class particle
{
public:
	float mass;
	float velocity;
	float xPosn;
	float yPosn;
	float force;
	particle()
	{
		mass = MASS;
		velocity = INITVELOCITY;
		float xPosn = XMIN;
		float yPosn = YMIN;
		force = 0.0;
	}

	void generateXPosn()
	{
		float random = ((float)rand()) / (float)RAND_MAX;
		float diff = XDIM - XMIN;
		float r = random * diff;
		xPosn = r;
	}

	void generateYPosn()
	{
		float random = ((float)rand()) / (float)RAND_MAX;
		float diff = YDIM - YMIN;
		float r = random * diff;
		yPosn = r;
	}
};

class node
{
public:
	node *nwChild;
	node *neChild;
	node *swChild;
	node *seChild;
	float xDim[2];
	float yDim[2];
	float centerOfMass;
	float cMassX;
	float cMassY;
	particle *content; //Will be set only and only if the external node
	bool isExternalNode;
	node *parentNode;
	vector<int> subParticleList;

	node()
	{
		nwChild = NULL;
		neChild = NULL;
		swChild = NULL;
		seChild = NULL;
		centerOfMass = 0.0;
		content = NULL;
		isExternalNode = false;
		parentNode = NULL;
		subParticleList.clear();
	}

	bool isLeafNode()
	{
		if (nwChild || neChild || swChild || seChild)
		{
			return false;
		}
		return true;
	}

	void getNWChild()
	{
		nwChild = new node();
		nwChild->xDim[0] = xDim[0];
		nwChild->xDim[1] = xDim[0] + ((xDim[1] - xDim[0]) / 2.0);
		nwChild->yDim[0] = yDim[0] + (yDim[1] - yDim[0]) / 2.0;
		nwChild->yDim[1] = yDim[1];
		nwChild->parentNode = this;
	}

	void getNEChild()
	{
		neChild = new node();
		neChild->xDim[0] = xDim[0] + ((xDim[1] - xDim[0]) / 2.0);
		neChild->xDim[1] = xDim[1];
		neChild->yDim[0] = yDim[0] + ((yDim[1] - yDim[0]) / 2.0);
		neChild->yDim[1] = yDim[1];
		neChild->parentNode = this;
	}

	void getSWChild()
	{
		swChild = new node();
		swChild->xDim[0] = xDim[0];
		swChild->xDim[1] = xDim[0] + ((xDim[1] - xDim[0]) / 2.0);
		swChild->yDim[0] = yDim[0];
		swChild->yDim[1] = yDim[0] + ((yDim[1] - yDim[0]) / 2.0);
		swChild->parentNode = this;
	}

	void getSEChild()
	{
		seChild = new node();
		seChild->xDim[0] = xDim[0] + ((xDim[1] - xDim[0]) / 2.0);
		seChild->xDim[1] = xDim[1];
		seChild->yDim[0] = yDim[0];
		seChild->yDim[1] = yDim[0] + ((yDim[1] - yDim[0]) / 2.0);
		seChild->parentNode = this;
	}
	void calcCenterOfMass() //Calculate the center of mass for the current node;
	{
	}
};


class barnesTree
{
public:
	node *m_root;
	long m_noOfParticles;
	long m_noOfSteps;
	float m_theta;
	particle **particleList;

	barnesTree(long noOfParticles, long noOfSteps, float theta)
	{
		m_noOfParticles = noOfParticles;
		m_noOfSteps = noOfParticles;
		m_theta = theta;
		m_root = new node();
		m_root->xDim[START] = XMIN;
		m_root->xDim[END] = XDIM;
		m_root->yDim[START] = YMIN;
		m_root->yDim[END] = XDIM;
		particleList = new particle*[m_noOfParticles];
	}

	//Member functions

	void prepareInitList()
	{
		particleList = new particle*[m_noOfParticles];
		srand(time(NULL));
		for (int i = 0; i < m_noOfParticles; i++)
		{
			particleList[i] = new particle();
			particleList[i]->generateXPosn();
			particleList[i]->generateYPosn();
			cout << "\n" << i << " xPosn = " << particleList[i]->xPosn << " yPosn = " << particleList[i]->yPosn;
		}
	}

	//Main engine of the Barnes Hut algorithm implementation
	bool run()
	{
		for (int i = 0; i < m_noOfParticles; i++)
		{
			insert(m_root,i);
		}

		calcCenterOfMass(m_root);

		for (int idx = 0; idx < m_noOfParticles; idx++)
		{
			calcForce(m_root,idx);
		}
		

		return true;
	}


	void insert(node *currNode,const int idx)
	{
		if (!currNode)
		{
			cout << "\n failed during insertion\n";
				return;
		}

		if (currNode->isLeafNode())
		{
			if (!currNode->content) //Place the particle in this node
			{
				currNode->content = particleList[idx];
				return;
			}
			else if (currNode->content)
			{
				splitNode(currNode, idx);
				currNode->content = NULL;
				insert(currNode, idx);
				return;
			}
		}
		else
		{
			//Check the quadrant to be moved into
			int region = checkQuadrant(currNode, particleList[idx]);
			switch (region)
			{
			case NW:
				if (currNode->nwChild)
				{
					insert(currNode->nwChild, idx);
				}
				else
				{
					//Set it as the nwChild of the currNode
					currNode->getNWChild();
					currNode->nwChild->content = particleList[idx];
				}
				break;
			case NE:
				if (currNode->neChild)
				{
					insert(currNode->neChild, idx);
				}
				else
				{
					//Set it as the neChild of the currNode
					currNode->getNEChild();
					currNode->neChild->content = particleList[idx];
				}
				break;
			case SW:
				if (currNode->swChild)
				{
					insert(currNode->swChild, idx);
				}
				else
				{
					//Set it as the nwChild of the currNode
					currNode->getSWChild();
					currNode->swChild->content = particleList[idx];
				}
				break;
			case SE:
				if (currNode->seChild)
				{
					insert(currNode->seChild, idx);
				}
				else
				{
					//Set it as the nwChild of the currNode
					currNode->getSEChild();
					currNode->seChild->content = particleList[idx];
				}
				break;
			default:
				cout << "\n Entered into the default option: Something went wrong";
				return;
			}
		}
	}

	void calcCenterOfMass(node *currNode)
	{
		if (currNode->isLeafNode())
		{
			currNode->centerOfMass = MASS;
			currNode->cMassX = currNode->content->xPosn;
			currNode->cMassY = currNode->content->yPosn;
			return;
		}
		float totalMass = 0.0;
		float cxPosn = 0.0;
		float cyPosn = 0.0;
		if (currNode->nwChild)
		{
			calcCenterOfMass(currNode->nwChild);
			totalMass += currNode->nwChild->centerOfMass;
			cxPosn += (currNode->nwChild->centerOfMass*currNode->nwChild->cMassX);
			cyPosn += (currNode->nwChild->centerOfMass*currNode->nwChild->cMassY);
		}

		if (currNode->neChild)
		{
			calcCenterOfMass(currNode->neChild);
			totalMass += currNode->neChild->centerOfMass;
			cxPosn += (currNode->neChild->centerOfMass*currNode->neChild->cMassX);
			cyPosn += (currNode->neChild->centerOfMass*currNode->neChild->cMassY);

		}

		if (currNode->swChild)
		{
			calcCenterOfMass(currNode->swChild);
			totalMass += currNode->swChild->centerOfMass;
			cxPosn += (currNode->swChild->centerOfMass*currNode->swChild->cMassX);
			cyPosn += (currNode->swChild->centerOfMass*currNode->swChild->cMassY);
		}

		if (currNode->seChild)
		{
			calcCenterOfMass(currNode->seChild);
			totalMass += currNode->seChild->centerOfMass;
			cxPosn += (currNode->seChild->centerOfMass*currNode->seChild->cMassX);
			cyPosn += (currNode->seChild->centerOfMass*currNode->seChild->cMassY);
		}
		currNode->centerOfMass = totalMass;
		currNode->cMassX = cxPosn / totalMass;
		currNode->cMassY = cyPosn / totalMass;
	}

	void calcForce(node *currNode,long int idx)
	{
		//(s/d<Theta)->Where where s is the width of the region represented by the internal node, and d is the 
		//distance between the body and the node's center-of-mass
		if (currNode->isLeafNode())
		{
			if (currNode->content && (currNode->content != particleList[idx])
			{
				particleList[idx]->force +=
			}
		}

	}
	int checkQuadrant(node *currNode, particle *temp)
	{
		//Check the NW quadrant of the current node
		float xMin = currNode->xDim[0];
		float xMax = currNode->xDim[0] + ((currNode->xDim[1] - currNode->xDim[0]) / 2.0);
		float yMin = currNode->yDim[0] + (currNode->yDim[1] - currNode->yDim[0]) / 2.0;
		float yMax = currNode->yDim[1];

		if (temp->xPosn >= xMin && temp->xPosn <= xMax)
		{
			if (temp->yPosn > yMin && temp->yPosn <= yMax)
			{
				return NW;
			}
		}

		//Check the NE Quadrant of the current node
		
		xMin = currNode->xDim[0] + ((currNode->xDim[1] - currNode->xDim[0]) / 2.0);
		xMax = currNode->xDim[1];
		yMin = currNode->yDim[0] + ((currNode->yDim[1] - currNode->yDim[0]) / 2.0);
		yMax = currNode->yDim[1];

		if (temp->xPosn > xMin && temp->xPosn <= xMax)
		{
			if (temp->yPosn > yMin && temp->yPosn <= yMax)
			{
				return NE;
			}
		}

		//Check the SW Quadrant of the current node
		xMin = currNode->xDim[0];
		xMax = currNode->xDim[0] + ((currNode->xDim[1] - currNode->xDim[0]) / 2.0);
		yMin = currNode->yDim[0];
		yMax = currNode->yDim[0] + ((currNode->yDim[1] - currNode->yDim[0]) / 2.0);

		if (temp->xPosn >= xMin && temp->xPosn <= xMax)
		{
			if (temp->yPosn >= yMin && temp->yPosn <= yMax)
			{
				return SW;
			}
		}

		//Check the SE Quadrant of the current node
		xMin = currNode->xDim[0] + ((currNode->xDim[1] - currNode->xDim[0]) / 2.0);
		xMax = currNode->xDim[1];
		yMin = currNode->yDim[0];
		yMax = currNode->yDim[0] + ((currNode->yDim[1] - currNode->yDim[0]) / 2.0);

		if (temp->xPosn > xMin && temp->xPosn <= xMax)
		{
			if (temp->yPosn >= yMin && temp->yPosn <= yMax)
			{
				return SE;
			}
		}
		//Must not reach till this point
		//If here that means an error
		return INVALID;
	}

	void splitNode(node *currNode, int idx)
	{
		//Here we can assume that the currNode is the leaf node but already has the content
		//So we need to move both the contents down
		if (!currNode)
		{
			cout << "\n Error::splitNode() got currNode as null-> returning without doing anythig";
			return;
		}

		if (currNode->xDim[1] - currNode->xDim[0] == 1.0f)
		{
			//Cannot split this node further.
			cout << "\n Cannot split this node further.";
			return;
		}
		int contentQuad = checkQuadrant(currNode, currNode->content);
		switch (contentQuad)
		{
		case NW:
			currNode->getNWChild();
			currNode->nwChild->content = currNode->content;
			break;
		case NE:
			currNode->getNEChild();
			currNode->neChild->content = currNode->content;
			break;
		case SW:
			currNode->getSWChild();
			currNode->swChild->content = currNode->content;
			break;
		case SE:
			currNode->getSEChild();
			currNode->seChild->content = currNode->content;
			break;
		}
	}
};

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		cout << "\n wrong usage";
		exit(0);
	}

	int noOfParticle = atoi(argv[1]);
	srand(time(NULL));
	barnesTree *tree = new barnesTree(noOfParticle, 4, THETA);

	for (int i = 0; i < noOfParticle; i++)
	{
		tree->particleList[i] = new particle();
		tree->particleList[i]->generateXPosn();
		tree->particleList[i]->generateYPosn();
		cout << "\n" << i << " xPosn = " << tree->particleList[i]->xPosn << " yPosn = " << tree->particleList[i]->yPosn;
	}

	//Build the tree now;
	bool retFlag = tree->run();
	if (retFlag == false)
	{
		cout << "\n failed to run the barnes hut tree";
		exit(1);
	}
	_getch();
	return 0;
}