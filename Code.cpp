#include <iostream> 
#include <cmath>
#include <math.h> 
#include<fstream>
#include<sstream>
#include <string> 
#include <stdlib.h>     /* srand, rand */
#include <vector>
#include <random>
#include <chrono>
#include <tuple>
#include <queue>
#include <time.h>
#include <stack> 
#include <list>
#include<ilcplex/ilocplex.h>

ILOSTLBEGIN

# define M_PI           3.14159265358979323846  /* pi */

using namespace std;
using namespace std::chrono;

struct Struct
{
	int p;
	int b;
	int value;
	bool temp[210];
};

struct ArrayDP
{
	double p;
	double b;
	bool temp[210];
};

struct Users
{
	double x;
	double y;
	double id;
	double md;
	double w;
	double p;
	double b;
	double total; /* summation of p and b */
	double max_p;
	double max_b;

	int index; /* to specify the number of users*/
};

struct Cloudlets
{
	double x;
	double y;
	double p;
	double b;
	int index;
};

struct values
{
	double delay;
	int index;
};
double ComputationTime = 0;
double MigrationTime = 0;
double OffloadingTime = 0;



double epsilon = 0.9; /* epsilon as an input ?*/

int user_total = 210; /* determines number of users */
int cloudlet_total = 35; /* determines number of cloudlets */

bool temp1[210];  /* user_total */
bool temp2[210];  /* user_total */

bool temp3[210]; /* inside Assignment function*/
double** UsageP = new double*[cloudlet_total]; /* user_total is the number of users */

double** UsageB = new double*[cloudlet_total]; /* user_total is the number of users */

ArrayDP* Array = new ArrayDP[user_total*10*user_total];


double total_size = 0;
int w = 5; /* window size*/
double gamma = 0.3; /* gamma -- weight */
double theta = 0.005*0.67;  /*  propagation time (IN MICROSECOND) per each unit of distance (meter) */
int DPRAC_Limit = 10; /* for DPRAC */
int load_balance = 0; /* for load balancing */
Users* temporary = new Users[w];  /* window size -- used in function PREDICT */
Users** temporary_test = new Users*[user_total];  /* used in SpecPredict */

Users** first_prediction = new Users*[user_total];
Users** second_prediction = new Users*[user_total];
Users** third_prediction = new Users*[user_total];
Users** fourth_prediction = new Users*[user_total];
Users** fifth_prediction = new Users*[user_total];


int** V = new int*[user_total]; /* V is used in OAMC and changes in every t */
double** r = new double*[user_total]; /* r is used in OAMC and changes in every t -- Rounding */
int** Temp_V = new int*[user_total]; /* V is used in OAMC and changes in every t */
int** F = new int*[user_total]; /* F is used in Assignment */
int** N = new int*[user_total]; /* N is used in Assignment */

int* Migrations = new int[user_total]; /* records the number of migrations of each application */
double* Latency = new double[user_total]; /* records the latency of each application in its lifetime */

int migration = 0;
int migration1 = 0;

int MigrationNewYork = 0;
int MigrationOrlando = 0;
int MigrationDACT = 0;

double min_x = INT_MAX;
double min_y = INT_MAX;
double max_x = INT_MIN;
double max_y = INT_MIN;
int T = INT_MAX; /* is the intersection of number of time slots */
Users** users = new Users*[user_total]; /* user_total is the number of users */
Cloudlets* cloudlets = new Cloudlets[cloudlet_total];
Cloudlets* C = new Cloudlets[cloudlet_total]; /* used in OAMC */
int* count_time = new int[user_total]; /* shows how many time slots each user runs */
int** R = new int*[user_total]; /* holds the index of the assigned cloudlet */

bool** S = new bool*[cloudlet_total]; /* holds the index of the assigned cloudlet */
bool* I = new bool[user_total]; /* used in DPRAC */

Users** First_Predict = new Users*[user_total]; /* user_total is the number of users */
Users** Second_Predict = new Users*[user_total]; /* user_total is the number of users */
Users** Third_Predict = new Users*[user_total]; /* user_total is the number of users */
Users** Fourth_Predict = new Users*[user_total]; /* user_total is the number of users */
Users** Fifth_Predict = new Users*[user_total]; /* user_total is the number of users */

bool flag_test = false;

typedef IloArray<IloBoolVarArray> BoolVarMatrix;
typedef IloArray<BoolVarMatrix>   BoolVar3Matrix;
typedef IloArray<BoolVar3Matrix>   BoolVar4Matrix;

double distance(double lat1, double lon1, double lat2, double lon2)
{
	return sqrt(pow(lat2 - lat1, 2) +
		pow(lon2 - lon1, 2) * 1.0);
}
static double haversine(double lat1, double lon1, double lat2, double lon2)
{
	/* distance between latitudes and longitudes */
	double dLat = (lat2 - lat1) * M_PI / 180.0;
	double dLon = (lon2 - lon1) *  M_PI / 180.0;

	/* convert to radians */
	lat1 = (lat1)* M_PI / 180.0;
	lat2 = (lat2)* M_PI / 180.0;

	/* apply formulae */
	double a = pow(sin(dLat / 2), 2) + pow(sin(dLon / 2), 2) * cos(lat1) * cos(lat2);
	double rad = 6371;
	double c = 2 * asin(sqrt(a));
	return rad * c;
}

void Initialize()
{
	for (int i = 0;i < user_total;i++)
	{
		Migrations[i] = 0;
		Latency[i] = 0;
	}
}

void Test()
{

	int d = 0;
	for (int i = 0;i<user_total;i++)
		d += Migrations[i];

	double d1 = 0;
	for (int i = 0;i<user_total;i++)
		d1 += Latency[i];
	cout << "migration  " << d << "  latency  " << d1 << endl;
}
void Print()
{
	for (int i = 0;i < user_total;i++)
	{
		cout<<"migration "<<Migrations[i]<<"  Latency "<<Latency[i]<<endl;
	}
}
bool compareTwoApplicationsBFD(Users user_1, Users user_2)
{
	return user_1.total >= user_2.total;
}
bool compareTwoCloudletsBFD(Cloudlets cloudlet_1, Cloudlets cloudlet_2)
{
	return cloudlet_1.p >= cloudlet_2.p && cloudlet_1.b >= cloudlet_2.b; 
}
bool compareTwoDelays(values value_1, values value_2)
{
	return value_1.delay <= value_2.delay;
}
void Seperator(string str, int User_Num, int timeslot)
{

	stringstream s(str); /* Used for breaking words */
	string word; /* to store individual words */

	int count = 0;
	while (s >> word)
	{

		stringstream geek(word);

		double read = 0;
		geek >> read;

		cout.precision(word.size());

		if (count == 0) /* latitude */
		{
			users[User_Num][timeslot - 1].x = read;
			
				if (read <= min_x)
					min_x = read;
				if (read>max_x)
					max_x = read;
					
		}
		if (count == 1) /* longitude */
		{
			users[User_Num][timeslot - 1].y = read;
			
			if (read <= min_y)
				min_y = read;
			if (read>max_y)
				max_y = read;

		}

		count++;
	}

}
void reading()
{

	int Num_User = 0; /* determines the number of the user */

	while (Num_User<user_total) /* user_total is the total number of users */
	{
		string  filename;
		stringstream out;
		out << Num_User + 1;
		filename = "C:\\Users\\Erfan\\Desktop\\" + out.str() + ".txt";

		string STRING;
		ifstream infile;
		infile.open(filename.c_str());
		int timeslots = 0;
		while (!infile.eof()) // To get you all the lines.
		{
			getline(infile, STRING); // Saves the line in STRING. 
			timeslots++;
		}
		infile.close();


		users[Num_User] = new Users[timeslots];  /* initializing the user */
		count_time[Num_User] = timeslots; /* count the number of active timeslots for each user*/
		R[Num_User] = new int[timeslots]; /* holds the assigned cloudlet of the user at each time slot */


		ifstream infile1;
		infile1.open(filename.c_str());
		while (!infile1.eof()) // To get you all the lines.
		{
			getline(infile1, STRING); // Saves the line in STRING. 
			Seperator(STRING, Num_User, timeslots);
			timeslots--;
		}
		infile1.close();

		Num_User++;
	}

}
void UserDistribution()
{
	for (int i = 0;i<user_total;i++)
	{

		/* for predicting idi^t: finding average*/
		std::random_device rd1;
		std::mt19937 gen1(rd1());
		std::uniform_real_distribution<> dis1(1*8, 8*20); /* in MB */  
		double sample1 = dis1(gen1);

		/* for predicting idi^t: finding variance*/
		std::random_device rd2;
		std::mt19937 gen2(rd2());
		std::uniform_real_distribution<> dis2(0, 2*8); /* in MB */  
		double sample2 = dis2(gen2);

		/* for predicting mdi^t: finding average*/
		std::random_device rd3;
		std::mt19937 gen3(rd3());
		std::uniform_real_distribution<> dis3(1*8, 3*16); /* in MB */ 
		double sample3 = dis3(gen3);

		/* for predicting mdi^t: finding variance*/
		std::random_device rd4;
		std::mt19937 gen4(rd4());
		std::uniform_real_distribution<> dis4(0, 1*8); /* in MB */ 
		double sample4 = dis4(gen4);

		/* for predicting wi^t: finding average*/
		std::random_device rd5;
		std::mt19937 gen5(rd5());
		std::uniform_real_distribution<> dis5(1, 80); /* in number of cycles needed (IN MILLIONS)*/ 
		double sample5 = dis5(gen5);

		/* for predicting wi^t: finding variance*/
		std::random_device rd6;
		std::mt19937 gen6(rd6());
		std::uniform_real_distribution<> dis6(0, 20); /* in number of cycles needed */  
		double sample6 = dis6(gen6);

		/* for predicting pi^t: finding average*/
		std::random_device rd7;
		std::mt19937 gen7(rd7());
		std::uniform_real_distribution<> dis7(1400, 1400); /* in M number of cycles per second */ 
		double sample7 = dis7(gen7);

		/* for predicting pi^t: finding variance*/
		std::random_device rd8;
		std::mt19937 gen8(rd8());
		std::uniform_real_distribution<> dis8(4400, 4400); /* in M number of cycles per second */ 
		double sample8 = dis8(gen8);

		/* for predicting bi^t: finding average*/
		std::random_device rd9;
		std::mt19937 gen9(rd9());
		std::uniform_real_distribution<> dis9(750, 1500); /* in MBPS bandwidth */ 
		double sample9 = dis9(gen9);

		/* for predicting bi^t: finding variance*/
		std::random_device rd10;
		std::mt19937 gen10(rd10());
		std::uniform_real_distribution<> dis10(100, 200); /* in MBPS bandwidth */ 
		double sample10 = dis10(gen10);

		for (int t = 0; t < count_time[i]; ++t)
		{
			bool flag = true;
			while (flag)
			{
				std::random_device rd1;
				std::mt19937 gen1(rd1());
				std::normal_distribution<double> d1(sample1/*mean[i]*/, sample2/*variance[i]*/);

				double sample11 = d1(gen1);
				users[i][t].id = sample11;

				users[i][t].index = i; /* to specify the number of users -- used in BFD and OAMC_BFD */

				std::random_device rd2;
				std::mt19937 gen2(rd2());
				std::normal_distribution<double> d2(sample3/*mean[i]*/, sample4/*variance[i]*/);

				double sample12 = d2(gen2);
				users[i][t].md = sample12;

				std::random_device rd3;
				std::mt19937 gen3(rd3());
				std::normal_distribution<double> d3(sample5/*mean[i]*/, sample6/*variance[i]*/);

				double sample13 = d3(gen3);
				users[i][t].w = sample13;

				std::random_device rd4;
				std::mt19937 gen4(rd4());
				std::normal_distribution<double> d4(sample7/*mean[i]*/, sample8/*variance[i]*/);

				double sample14 = d4(gen4);
				users[i][t].p = sample14;

				std::random_device rd5;
				std::mt19937 gen5(rd5());
				std::normal_distribution<double> d5(sample9/*mean[i]*/, sample10/*variance[i]*/);

				double sample15 = d5(gen5);
				users[i][t].b = sample15;

				if (users[i][t].p != 0 && users[i][t].b != 0 && users[i][t].id>=0 && users[i][t].md>0 && users[i][t].w >= 0 && users[i][t].p>0 && users[i][t].b>0)
					flag = false;
				
			}

		}
	}


}
void CloudletDistribution()
{

	for (int i = 0; i < cloudlet_total; ++i)
	{
		/* for x of the clodulet:*/
		std::random_device rd1;
		std::mt19937 gen1(rd1());
		std::uniform_real_distribution<> dis1(min_x, max_x); /* latitude */
		double sample1 = dis1(gen1);

		/* for y of the clodulet:*/
		std::random_device rd2;
		std::mt19937 gen2(rd2());
		std::uniform_real_distribution<> dis2(min_y, max_y); /* latitude */
		double sample2 = dis2(gen2);

		/*for P U(1,20) */ 		/* for p of the clodulet:*/
		std::random_device rdP;
		std::mt19937 genP(rdP());
		std::uniform_real_distribution<> disP(1, 20); /* number of VMs cumulative P */
		int sampleP = disP(genP);

		int counter = 0;
		double sample3 = 0;
		while (counter < sampleP)
		{

		   std::random_device rd3;
		   std::mt19937 gen3(rd3());
		   std::normal_distribution<double> dis3(25000,50000);
			
		   double test = dis3(gen3);
		   if(test>0)
			 {sample3 += test; counter++;}

		}
		

		/* for b of the clodulet:*/
		std::random_device rd4;
		std::mt19937 gen4(rd4());
		std::uniform_real_distribution<> dis4(5000, 10000); /* b */  
		double sample4 = dis4(gen4);

		cloudlets[i].x = sample1;
		cloudlets[i].y = sample2;
		cloudlets[i].p = sample3;
		cloudlets[i].b = sample4;


	}


}
void Dataset() /* set up a dataset */
{
	reading(); /* for reading lat and lon of the users */
	UserDistribution(); /* generating other user specifications in their lifetime such as b,p,w,id,md using distribution*/
	CloudletDistribution(); /* generating cloudlet specifications such as x,y,p,b */
}

Users* PREDICT_NEW(int i, int t)
{
	temporary[0].x = First_Predict[i][t].x;
	temporary[0].id = First_Predict[i][t].id;
	temporary[0].y = First_Predict[i][t].y;
	temporary[0].md = First_Predict[i][t].md;
	temporary[0].w = First_Predict[i][t].w;
	temporary[0].p = First_Predict[i][t].p;
	temporary[0].b = First_Predict[i][t].b;

	temporary[1].x = Second_Predict[i][t].x;
	temporary[1].id = Second_Predict[i][t].id;
	temporary[1].y = Second_Predict[i][t].y;
	temporary[1].md = Second_Predict[i][t].md;
	temporary[1].w = Second_Predict[i][t].w;
	temporary[1].p = Second_Predict[i][t].p;
	temporary[1].b = Second_Predict[i][t].b;

	temporary[2].x = Third_Predict[i][t].x;
	temporary[2].id = Third_Predict[i][t].id;
	temporary[2].y = Third_Predict[i][t].y;
	temporary[2].md = Third_Predict[i][t].md;
	temporary[2].w = Third_Predict[i][t].w;
	temporary[2].p = Third_Predict[i][t].p;
	temporary[2].b = Third_Predict[i][t].b;

	temporary[3].x = Fourth_Predict[i][t].x;
	temporary[3].id = Fourth_Predict[i][t].id;
	temporary[3].y = Fourth_Predict[i][t].y;
	temporary[3].md = Fourth_Predict[i][t].md;
	temporary[3].w = Fourth_Predict[i][t].w;
	temporary[3].p = Fourth_Predict[i][t].p;
	temporary[3].b = Fourth_Predict[i][t].b;

	temporary[4].x = Fifth_Predict[i][t].x;
	temporary[4].id = Fifth_Predict[i][t].id;
	temporary[4].y = Fifth_Predict[i][t].y;
	temporary[4].md = Fifth_Predict[i][t].md;
	temporary[4].w = Fifth_Predict[i][t].w;
	temporary[4].p = Fifth_Predict[i][t].p;
	temporary[4].b = Fifth_Predict[i][t].b;

	return temporary;
}

Users* PREDICT(int i, int t)
{
	for (int h = 0;h < w;h++)
	{
		/* for predicting x */
		temporary[h].x = users[i][t].x;
		double dif = users[i][t].x - users[i][1].x;
		if (t != 0)
			temporary[h].x += (double(h + 1) / double(t))*dif;

		/* for predicting id */
		double sum1 = 0;
		for (int counter = 1; counter <= t;counter++)
			sum1 += (counter)*users[i][counter].id;
		sum1 /= (t )*(t + 1) / 2;
		temporary[h].id = sum1;

		/* for predicting y */
		temporary[h].y = users[i][t].y;
		double dif1 = users[i][t].y - users[i][1].y;
		if (t != 0)
			temporary[h].y += (double(h + 1) / double(t))*dif1;

		/* for predicting md */
		int sum2 = 0;
		for (int counter = 1; counter <= t;counter++)
			sum2 += (counter)*users[i][counter].md;
		sum2 /= ((t )*(t + 1) / 2);
		temporary[h].md = sum2;

		/* for predicting w */
		int sum3 = 0;
		for (int counter = 1; counter <= t;counter++)
			sum3 += (counter)*users[i][counter].w;
		sum3 /= ((t)*(t + 1) / 2);
		temporary[h].w = sum3;

		/* for predicting p */
		int sum4 = 0;
		for (int counter = 1; counter <= t;counter++)
			sum4 += (counter)*users[i][counter].p;
		sum4 /= ((t)*(t + 1) / 2);
		temporary[h].p = sum4;

		/* for predicting b */
		int sum5 = 0;
		for (int counter = 1; counter <= t;counter++)
			sum5 += (counter)*users[i][counter].b;
		sum5 /= ((t)*(t + 1) / 2);
		temporary[h].b = sum5;

	}
	return temporary;
}

 void DPRAC_NEW(int** V, Cloudlets* C1, int j, int t)
{
	/* initializing I */

	for (int i = 0;i<user_total;i++)
		I[i] = false;

	/*make the values positive */
	int min = INT_MAX;
	for (int i = 0;i<user_total;i++)
		if (V[i][j]!=INT_MAX && V[i][j] <= min)
			min = V[i][j];
	if (min<0)
		for (int i = 0;i<user_total;i++)
			if (V[i][j] != INT_MAX)
			V[i][j] += (-min + 1);

	int M = INT_MIN;
	for (int i = 0;i<user_total;i++)
		if (V[i][j] >= M && V[i][j]!=INT_MAX && V[i][j] >= 0)
			M = V[i][j];
	int Mi = INT_MAX;
	for (int i = 0;i<user_total;i++)
		if (V[i][j]!=INT_MAX && V[i][j] <= Mi && V[i][j] >= 0)
			Mi = V[i][j];


    /* make the values starting from 0*/
	for (int i = 0;i<user_total;i++)
	  if(V[i][j] != INT_MAX && V[i][j]>=0)
		V[i][j] = M+1-V[i][j];

	int M2 = INT_MIN;
	for (int i = 0;i<user_total;i++)
		if (V[i][j] != INT_MAX && V[i][j] >= M2 && V[i][j] >= 0)
			M2 = V[i][j];
	int Mi2 = INT_MAX;
	for (int i = 0;i<user_total;i++)
		if (V[i][j] != INT_MAX && V[i][j] <= Mi2 && V[i][j] >= 0)
			Mi2 = V[i][j];
	

	//double epsilon = 0.5; /* epsilon should be an input */

	double lambda = (epsilon*M2) / user_total;
	if(lambda < 1)
	 lambda = 1;


	if (lambda < 0)
	{
		cout << "negative " << lambda << endl;
	}


	for(int i=0;i<user_total;i++)
	  for(int j=0;j<cloudlet_total;j++)
		  if (V[i][j] != INT_MAX && V[i][j] >= 0)
	        r[i][j] = (double)V[i][j]/ lambda;
		  else 
		    r[i][j] = INT_MAX;

	double rmax = INT_MIN;
	for (int i = 0;i<user_total;i++)
		if (r[i][j]  >= rmax && V[i][j]!=INT_MAX && V[i][j] >= 0)
			rmax = r[i][j];

	double rmin = INT_MAX;
	for (int i = 0;i<user_total;i++)
		if (r[i][j] <= rmin && V[i][j]!=INT_MAX && V[i][j] >= 0)
			rmin = r[i][j];


	int tst = 0;


	int size = M2/ lambda; /* Size is v* */


	int count = 0;
	for(int i=0;i<user_total;i++)
		if(R[i][t]==-1 && V[i][j] != INT_MAX && V[i][j] >= 0)
		   count++;
		
	int * indices = new int[count];
	

	int count2 = 0;
	for (int i = 0;i<user_total;i++)
		if (R[i][t] == -1 && V[i][j] != INT_MAX && V[i][j] >= 0)
			{indices[count2]=i;count2++;}

	int ArraySize = size*count;
		

	for (int i = 1;i < count*size;i++)
	{
		if (size!=0)
		{
			Array[i].p = 1000000;
			Array[i].b = 1000000;
			for (int j = 0;j<user_total;j++)
				Array[i].temp[j] = 0;
		}
			
	}

	Array[0].p = 0;
	Array[0].b = 0;
	for (int j = 0;j<user_total;j++)
		Array[0].temp[j] = 0;


	for(int i=1;i<=count;i++)
		for (int v = (i - 1)*size;v >= 0;v--)
		{
			

				int value = V[indices[i-1]][j]/ lambda;

				
				 if (Array[v + value].p > (Array[v].p + users[indices[i-1]][t].p))
				 {
				 
					 Array[v + value].p = Array[v].p + users[indices[i-1]][t].p;
										
					 for (int j = 0;j<user_total;j++)
						 Array[v + value].temp[j] = Array[v].temp[j];

					 Array[v + value].temp[indices[i-1]] = 1;
				 }
				 if ( Array[v + value].b > (Array[v].b + users[indices[i-1]][t].b))
				 {

					 Array[v + value].b = Array[v].b + users[indices[i-1]][t].b;

					 for (int j = 0;j<user_total;j++)
						 Array[v + value].temp[j] = Array[v].temp[j];

					 Array[v + value].temp[indices[i - 1]] = 1;
				 }

		} 
	
        double max = INT_MIN;
		int index = -1;
		for (int v = count*size-1;v >=0; v--)
		{
			

			 if (Array[v].p <= C1[j].p && Array[v].b <= C1[j].b && Array[v].p>0 && Array[v].b>0)
			 {
				index = v;
				v = -1;
			 }
			
		}

		if (count == 1 && index==-1)
		{
			int value = V[indices[0]][j] / lambda;
			Array[value].p = users[indices[0]][t].p;

			Array[value].b = users[indices[0]][t].b;

			for (int j = 0;j<user_total;j++)
				Array[value].temp[j] = 0;

			Array[value].temp[indices[0]] = 1;
		   if (Array[value].p <= C1[j].p && Array[value].b <= C1[j].b && Array[value].p>0 && Array[value].b>0)
			  index = value;
		}

		for (int i = 0;i < user_total;i++)
		{
		    if(index==-1)
			  temp3[i] = 0;
			else
			  temp3[i] = Array[index].temp[i];
		}
			
}

void DPRAC_Sample(int** V, Cloudlets* C1, int j, int t)
{
	/* initializing I */

	for (int i = 0;i<user_total;i++)
		I[i] = false;


	/*make the values positive */
	int min = INT_MAX;
	for (int i = 0;i<user_total;i++)
		if (V[i][j] != INT_MAX && V[i][j] <= min)
			min = V[i][j];
	if (min<0)
		for (int i = 0;i<user_total;i++)
			if (V[i][j] != INT_MAX)
				V[i][j] += (-min + 1);

	values* temporary = new values[user_total];
	for (int i = 0;i<user_total;i++)
	{
		temporary[i].delay = V[i][j];
		temporary[i].index = i;
	}
	std::sort(temporary, temporary + user_total, compareTwoDelays);

	int sample_check = 0;
	int M = INT_MIN;
	int Mi = INT_MAX;
	for (int i = 0;i<user_total;i++)
		if (V[temporary[i].index][j] != INT_MAX && V[temporary[i].index][j] >= 0  && R[temporary[i].index][t]==-1)
			{
			  sample_check++;
			  if (sample_check <= DPRAC_Limit)
			  {
				  if (V[temporary[i].index][j] >= M)
					  M = V[temporary[i].index][j];

				  if(V[temporary[i].index][j] <= Mi)
					  Mi = V[temporary[i].index][j];

			  }
			  
			}


	/* make the values starting from 0*/
	int sample_check2 = 0;
	for (int i = 0;i<user_total;i++)
		if (V[temporary[i].index][j] != INT_MAX && V[temporary[i].index][j] >= 0 && R[temporary[i].index][t] == -1)
		{
			sample_check2++;
			if(sample_check2<=DPRAC_Limit)
			   V[temporary[i].index][j] = M + 1 - V[temporary[i].index][j];
		}

	
	int M2 = INT_MIN;
	int Mi2 = INT_MAX;
	int sample_check3 = 0;

	for (int i = 0;i<user_total;i++)
		if ( V[temporary[i].index][j] != INT_MAX && V[temporary[i].index][j] >= 0  && R[temporary[i].index][t]==-1)
			{
			  sample_check3++;
			  if (sample_check3 <= DPRAC_Limit)
			  {
				  if (V[temporary[i].index][j] >= M2)
					  M2 = V[temporary[i].index][j];

				  if(V[temporary[i].index][j] <= Mi2)
					  Mi2 = V[temporary[i].index][j];

			  }
			  
			}


	int user_counter = 0;
	for (int i = 0;i<user_total;i++)
		if (V[temporary[i].index][j] != INT_MAX && V[temporary[i].index][j] >= 0 && R[temporary[i].index][t]==-1 && user_counter<DPRAC_Limit)
		   user_counter++;


	double lambda = (epsilon*M2) / user_counter;
	if (lambda < 1)
		lambda = 1;


	if (lambda < 0)
	{
		cout << "negative " << lambda << endl;
	}

	for (int i = 0;i<user_total;i++)
		for (int j = 0;j<cloudlet_total;j++)
			if (V[temporary[i].index][j] != INT_MAX && V[temporary[i].index][j] >= 0)
				r[temporary[i].index][j] = (double)V[temporary[i].index][j] / lambda;
			else
				r[temporary[i].index][j] = INT_MAX;

	int sample_check4 = 0;
	double rmax = INT_MIN;
	double rmin = INT_MAX;

	 for (int i = 0;i<user_total;i++)
		if (V[temporary[i].index][j] != INT_MAX && V[temporary[i].index][j] >= 0  && R[temporary[i].index][t] == -1)
		{
			sample_check4++;
			if (sample_check4 <= DPRAC_Limit)
			{
				if (r[temporary[i].index][j] >= rmax)
					rmax = r[temporary[i].index][j];

				if (r[temporary[i].index][j] <= rmin) 
					rmin = r[temporary[i].index][j];

			}

		}



	int size = M2 / lambda; /* Size is v* */

	int count = 0;
	for (int i = 0;i<user_total;i++)
		if (R[temporary[i].index][t] == -1 && V[temporary[i].index][j] != INT_MAX && V[temporary[i].index][j] >= 0 && count<DPRAC_Limit)
			count++;
	int * indices = new int[count];
	
	int count2 = 0;
	for (int i = 0;i<user_total;i++)
		if (R[temporary[i].index][t] == -1 && V[temporary[i].index][j] != INT_MAX && V[temporary[i].index][j] >= 0 && count2<DPRAC_Limit)
		{
			indices[count2] = temporary[i].index;
			count2++;
		}


	int ArraySize = size*count;

	for (int i = 1;i < count*size;i++)
	{
		if (size != 0)
		{
			Array[i].p = 1000000;
			Array[i].b = 1000000;
			for (int j = 0;j<user_total;j++)
				Array[i].temp[j] = 0;
		}

	}
	
	Array[0].p = 0;
	Array[0].b = 0;
	for (int j = 0;j<user_total;j++)
		Array[0].temp[j] = 0;


	for (int i = 1;i <= count;i++)
		for (int v = (i - 1)*size;v >= 0;v--)
		{

			int value = V[indices[i - 1]][j] / lambda;


			if (Array[v + value].p > (Array[v].p + users[indices[i - 1]][t].p))
			{
				
				Array[v + value].p = Array[v].p + users[indices[i - 1]][t].p;

				for (int j = 0;j<user_total;j++)
					Array[v + value].temp[j] = Array[v].temp[j];

				Array[v + value].temp[indices[i - 1]] = 1;
			}
			if (Array[v + value].b >(Array[v].b + users[indices[i - 1]][t].b))
			{

				Array[v + value].b = Array[v].b + users[indices[i - 1]][t].b;

				for (int j = 0;j<user_total;j++)
					Array[v + value].temp[j] = Array[v].temp[j];

				Array[v + value].temp[indices[i - 1]] = 1;
			}

		}


	double max = INT_MIN;
	int index = -1;
	for (int v = count*size - 1;v >= 0; v--)
	{

		if ( Array[v].p <= C1[j].p && Array[v].b <= C1[j].b && Array[v].p>0 && Array[v].b>0)
		{
			index = v;
			v = -1;
		}

	}

	if (count == 1 && index == -1)
	{
		int value = V[indices[0]][j] / lambda;
		Array[value].p = users[indices[0]][t].p;

		Array[value].b = users[indices[0]][t].b;

		for (int j = 0;j<user_total;j++)
			Array[value].temp[j] = 0;

		Array[value].temp[indices[0]] = 1;
		if (Array[value].p <= C1[j].p && Array[value].b <= C1[j].b && Array[value].p>0 && Array[value].b>0)
			index = value;
	}

	for (int i = 0;i < user_total;i++)
	{
		if (index == -1)
			temp3[i] = 0;
		else
			temp3[i] = Array[index].temp[i];
	}

}


bool** Assignment(int** V, Cloudlets* C, int t)
{
	int j = 0;
	for (int i = 0;i<user_total;i++)
		for (int j = 0;j<cloudlet_total;j++)
			Temp_V[i][j] = V[i][j];

	while (j < cloudlet_total)
	{
		
		DPRAC_Sample(Temp_V, C, j, t);
		for (int num = 0;num<user_total;num++)
			S[j][num] = temp3[num];

		for (int i = 0;i < user_total;i++)
		{
			for (int k = 0;k < cloudlet_total;k++)
			{
				if (S[j][i] == 1 || k == j)
				{
					F[i][k] = Temp_V[i][j];
				}
				else
					F[i][k] = 0;

			}
		}


		for (int i = 0;i < user_total;i++)
		{
			for (int k = 0;k < cloudlet_total;k++)
			{
				N[i][k] = Temp_V[i][k] - F[i][k];
			}
		}

		for (int i = 0;i < user_total;i++)
		{
			for (int k = 0;k < cloudlet_total;k++)
			{
				Temp_V[i][k] = N[i][k]; /* without column of cj*/
										
			}
		}

		j = j + 1;

	}

	for (int i = cloudlet_total - 2;i >= 0;i--) /* lines 13 to 15 of the code */
	{
		//S[i][] = 
		for (int j = 0;j < user_total;j++)
		{
			if (S[i][j] == 1)
			{
				for (int k = i + 1;k<cloudlet_total;k++)
					if (S[k][j] == 1)
						S[i][j] = 0;
			}
		}


	}

	return S;
}

bool** AssignmentNEW(int** V, Cloudlets* C, int t)
{
	int j = 0;

	for (int i = 0;i<user_total;i++)
		for (int j = 0;j<cloudlet_total;j++)
			Temp_V[i][j] = V[i][j];

	while (j < cloudlet_total)
	{

		DPRAC_NEW(Temp_V, C, j, t);
		for (int num = 0;num<user_total;num++)
			S[j][num] = temp3[num];

		for (int i = 0;i < user_total;i++)
		{
			for (int k = 0;k < cloudlet_total;k++)
			{
				if (S[j][i] == 1 || k == j)
				{
					F[i][k] = Temp_V[i][j];
				}
				else
					F[i][k] = 0;

			}
		}


		for (int i = 0;i < user_total;i++)
		{
			for (int k = 0;k < cloudlet_total;k++)
			{
				N[i][k] = Temp_V[i][k] - F[i][k];
			}
		}


		for (int i = 0;i < user_total;i++)
		{
			for (int k = 0;k < cloudlet_total;k++)
			{
				Temp_V[i][k] = N[i][k]; /* without column of cj*/
										
			}
		}

		j = j + 1;

	}

	for (int i = cloudlet_total - 2;i >= 0;i--) /* lines 13 to 15 of the code */
	{
		for (int j = 0;j < user_total;j++)
		{
			if (S[i][j] == 1)
			{
				for (int k = i + 1;k<cloudlet_total;k++)
					if (S[k][j] == 1)
						S[i][j] = 0;
			}
		}


	}

	return S;
}


double OAMC()
{
	double D = 0; /* measures delay */

	for (int k = 0;k<user_total;k++)
		if (count_time[k] <= T)
			T = count_time[k];
	T = 20; 
	for (int t = 0;t <T;t++)
	{
		for (int i = 0;i < user_total;i++)
		{
			R[i][t] = -1; /* initialization*/

			Users* u = PREDICT_NEW(i, t); /* gets w predicted specifications including u[0], u[1], ..., u[w-1]*/

			for (int j = 0;j < cloudlet_total;j++) 	/* calculating v_{ij}^t */
			{
				V[i][j] = 0.001*(haversine(users[i][t].x, users[i][t].y, cloudlets[j].x, cloudlets[j].y) * 1000 * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000);
				
				for (int counter = 0;counter < w;counter++)
				{
					V[i][j] += 0.001*(pow(gamma, counter + 1)*(haversine(u[counter].x, u[counter].y, cloudlets[j].x, cloudlets[j].y) * 1000 * theta + (u[counter].id / u[counter].b + u[counter].w / u[counter].p)) * 1000000);
					
				}


				/* checking migration*/
				if (t >= 1 && j != R[i][t - 1])
				{
					int a = R[i][t - 1];
					V[i][j] += 0.001*(haversine(cloudlets[a].x, cloudlets[a].y, cloudlets[j].x, cloudlets[j].y) * 1000 * theta + (users[i][t].md / users[i][t].b) * 1000000);

				}

				if(V[i][j]<0) 
				   V[i][j] = INT_MAX;
			}

		}

		for (int i = 0;i<user_total; i++)
			for (int j = 0;j < cloudlet_total;j++)
		        V[i][j] = V[i][j]/500;

		for (int num = 0;num<cloudlet_total;num++)
		{
			C[num].b = cloudlets[num].b;
			C[num].p = cloudlets[num].p;
			C[num].x = cloudlets[num].x;
			C[num].y = cloudlets[num].y;
		}
		bool flag = true;
		while (flag)
		{
			
			Assignment(V, C, t); 

								

			/* mapping the result of Assignment into R */
			for (int a = 0;a<cloudlet_total;a++)
				for (int b = 0;b<user_total;b++)
					if (S[a][b] == 1)
					{
						R[b][t] = a;
		
						C[a].b -= users[b][t].b;
						C[a].p -= users[b][t].p;
				
					}
			flag = false;
			int count_unassigned = 0;
			for (int i = 0;i<user_total;i++)
				if (R[i][t] == -1) /* if there is an unassigned user */
				{
					flag = true;
					count_unassigned++;

				}

		}

		for (int i = 0;i < user_total;i++)
		{
			int index = R[i][t];
			D += haversine(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * 1000 * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;
		    Latency[i]+= haversine(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * 1000 * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;
		}

		if (t >= 1)
			for (int i = 0;i < user_total;i++)
			{

				int index1 = R[i][t - 1];
				int index2 = R[i][t];

				if (index1 != index2)
				{
				    Migrations[i]++;
					migration++; D += haversine(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * 1000 * theta + (users[i][t].md / users[i][t].b) * 1000000;
					Latency[i]+= haversine(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * 1000 * theta + (users[i][t].md / users[i][t].b) * 1000000;
				}
			}
	}

	return D;
}


double S-OAMC-WP()
{
	double D = 0; /* measures delay */

				 
	for (int k = 0;k<user_total;k++)
		if (count_time[k] <= T)
			T = count_time[k];
	T = 20; 
	for (int t = 1;t <=T;t++)
	{
		
		for (int i = 0;i < user_total;i++)
		{
			R[i][t] = -1; /* initialization*/
						 
			for (int j = 0;j < cloudlet_total;j++) 	/* calculating v_{ij}^t */
			{
				V[i][j] = 0.001*(distance(users[i][t].x, users[i][t].y, cloudlets[j].x, cloudlets[j].y) * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000);

				/* checking migration*/
				if (t >= 2 && j != R[i][t - 1])
				{
					int a = R[i][t - 1];
					V[i][j] += 0.001*(distance(cloudlets[a].x, cloudlets[a].y, cloudlets[j].x, cloudlets[j].y) * theta + (users[i][t].md / users[i][t].b) * 1000000);
					
				}

				if (V[i][j]<0)
					V[i][j] = INT_MAX;
			}

		}
		for (int i = 0;i<user_total; i++)
			for (int j = 0;j < cloudlet_total;j++)
				V[i][j] = V[i][j] / 500;

		for (int num = 0;num<cloudlet_total;num++)
		{
			C[num].b = cloudlets[num].b;
			C[num].p = cloudlets[num].p;
			C[num].x = cloudlets[num].x;
			C[num].y = cloudlets[num].y;
		}
		bool flag = true;
		while (flag)
		{

			Assignment(V, C, t); 

			/* mapping the result of Assignment into R */
			for (int a = 0;a<cloudlet_total;a++)
				for (int b = 0;b<user_total;b++)
					if (S[a][b] == 1)
					{
						R[b][t] = a;

						C[a].b -= users[b][t].b;
						C[a].p -= users[b][t].p;

					}
			flag = false;
			for (int i = 0;i<user_total;i++)
				if (R[i][t] == -1) /* if there is an unassigned user */
				{
					flag = true;
				}

		}

		/* TEST */
		bool* testarray = new bool[cloudlet_total];
		for (int i = 0;i<cloudlet_total;i++)
			testarray[i] = 0;

		for (int i = 0;i<user_total;i++)
		{
			int index;
			index = R[i][t];
			testarray[index] = 1;

			UsageP[index][t] += users[i][t].p;
			UsageB[index][t] += users[i][t].b;

		}
		int count = 0;
		for (int i = 0;i<cloudlet_total;i++)
			if (testarray[i] == 1)
				count++;

		load_balance += count;

		for (int i = 0;i < user_total;i++)
		{
			int index = R[i][t];
			D += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;
			Latency[i] += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;
			
			/* calculate metrics */
			OffloadingTime += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * theta + (users[i][t].id / users[i][t].b) * 1000000;
			ComputationTime += (users[i][t].w / users[i][t].p) * 1000000;
		}

		if (t >= 2)
			for (int i = 0;i < user_total;i++)
			{

				int index1 = R[i][t - 1];
				int index2 = R[i][t];

				if (index1 != index2)
				{
					if (i<user_total / 3)
						MigrationNewYork++;
					if (i >= user_total / 3 && i<(user_total * 2) / 3)
						MigrationOrlando++;
					if (i >= (user_total * 2) / 3 && i<user_total)
						MigrationDACT++;

				    Migrations[i]++;
					migration1++; D += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;
					Latency[i] += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;
					MigrationTime += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;
				}
			}
	}

	return D;
}

double S-OAMC()
{
	double D = 0; /* measures delay */

	T = 20; 

	for (int t = 1;t <=T;t++)
	{
		if (t%w == 1 || w==1)
		{


			for (int i = 0;i < user_total;i++)
			{
				R[i][t] = -1; /* initialization*/

				Users* u = PREDICT_NEW(i, t); /* gets w predicted specifications including u[0], u[1], ..., u[w-1]*/

				for (int j = 0;j < cloudlet_total;j++) 	/* calculating v_{ij}^t */
				{
				  
				    V[i][j] = 0.001*(distance(users[i][t].x, users[i][t].y, cloudlets[j].x, cloudlets[j].y) * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000);
									  
					int w1; /* w1 is minimum between w and 20-t*/
					if(20-t<w)
					  w1 = 20-t;
					else
					  w1 = w;
					for (int counter = 0;counter < w1;counter++)
					{
						V[i][j] += 0.001*(pow(gamma, counter + 1)*(distance(u[counter].x, u[counter].y, cloudlets[j].x, cloudlets[j].y) * theta + (u[counter].id / u[counter].b + u[counter].w / u[counter].p)) * 1000000);

					}


					/* checking migration*/
					if (t >= 2 && j != R[i][t - 1])
					{
						int a = R[i][t - 1];
						V[i][j] += 0.001*(distance(cloudlets[a].x, cloudlets[a].y, cloudlets[j].x, cloudlets[j].y) * theta + (users[i][t].md / users[i][t].b) * 1000000);

					}

					if (V[i][j]<0)
						V[i][j] = INT_MAX;
				  
				}

			}

			for (int i = 0;i<user_total; i++)
				for (int j = 0;j < cloudlet_total;j++)
					V[i][j] = V[i][j] / 500;
			
			for (int num = 0;num<cloudlet_total;num++)
			{
				C[num].b = cloudlets[num].b;
				C[num].p = cloudlets[num].p;
				C[num].x = cloudlets[num].x;
				C[num].y = cloudlets[num].y;
			}
			bool flag = true;
			while (flag)
			{
				
				Assignment(V, C, t); 
			    
				/* mapping the result of Assignment into R */
				for (int a = 0;a<cloudlet_total;a++)
					for (int b = 0;b<user_total;b++)
						if (S[a][b] == 1)
						{
							R[b][t] = a;

							C[a].b -= users[b][t].b;
							C[a].p -= users[b][t].p;
						
						}
				flag = false;
				for (int i = 0;i<user_total;i++)
					if (R[i][t] == -1) /* if there is an unassigned user */
					{
						flag = true;
					}

			}

		}

		else
		{
			for (int i = 0;i < user_total;i++)
			{
				R[i][t] = R[i][t - 1];
			}
		}

		/* TEST */
		bool* testarray = new bool[cloudlet_total];
		for (int i = 0;i<cloudlet_total;i++)
			testarray[i] = 0;

		for (int i = 0;i<user_total;i++)
		{
			int index;
			index = R[i][t];
			testarray[index] = 1;

			UsageP[index][t] += users[i][t].p;
			UsageB[index][t] += users[i][t].b;

		}
		int count = 0;
		for (int i = 0;i<cloudlet_total;i++)
			if (testarray[i] == 1)
				count++;

		load_balance += count;

		/* calculating delay */
		for (int i = 0;i < user_total;i++)
		{
			int index = R[i][t];
			D += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;
			Latency[i] += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;

			/* calculate metrics */
			OffloadingTime += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * theta + (users[i][t].id / users[i][t].b) * 1000000;
			ComputationTime += (users[i][t].w / users[i][t].p) * 1000000;


		}

		if (t >= 2)
			for (int i = 0;i < user_total;i++)
			{

				int index1 = R[i][t - 1];
				int index2 = R[i][t];

				if (index1 != index2)
				{
					if (i<user_total / 3)
						MigrationNewYork++;
					if (i >= user_total / 3 && i<(user_total * 2) / 3)
						MigrationOrlando++;
					if (i >= (user_total * 2) / 3 && i<user_total)
						MigrationDACT++;

				    Migrations[i]++;
					migration++;  D += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;
					Latency[i] += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;

					MigrationTime += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;
				}
			}
	}

	return D;
}

double OAMC()
{
	double D = 0; /* measures delay */

	T = 20; 

	for (int t = 1;t <= T;t++)
	{
		if (t%w == 1 || w==1)
		{


			for (int i = 0;i < user_total;i++)
			{
				R[i][t] = -1; /* initialization*/
							 
				Users* u = PREDICT_NEW(i, t); /* gets w predicted specifications including u[0], u[1], ..., u[w-1]*/
											  

				for (int j = 0;j < cloudlet_total;j++) 	/* calculating v_{ij}^t */
				{

					V[i][j] = 0.001*(distance(users[i][t].x, users[i][t].y, cloudlets[j].x, cloudlets[j].y) * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000);

					int w1; /* w1 is minimum between w and 20-t*/
					if (20 - t<w)
						w1 = 20 - t;
					else
						w1 = w;
					for (int counter = 0;counter < w1;counter++)
					{
						V[i][j] += 0.001*(pow(gamma, counter + 1)*(distance(u[counter].x, u[counter].y, cloudlets[j].x, cloudlets[j].y) * theta + (u[counter].id / u[counter].b + u[counter].w / u[counter].p)) * 1000000);
					
					}

					/* checking migration*/
					if (t >= 2 && j != R[i][t - 1])
					{
						int a = R[i][t - 1];
						V[i][j] += 0.001*(distance(cloudlets[a].x, cloudlets[a].y, cloudlets[j].x, cloudlets[j].y) * theta + (users[i][t].md / users[i][t].b) * 1000000);

					}

					if (V[i][j]<0)
						V[i][j] = INT_MAX;

				}

			}

			for (int i = 0;i<user_total; i++)
				for (int j = 0;j < cloudlet_total;j++)
					V[i][j] = V[i][j] / 500;

			for (int num = 0;num<cloudlet_total;num++)
			{
				C[num].b = cloudlets[num].b;
				C[num].p = cloudlets[num].p;
				C[num].x = cloudlets[num].x;
				C[num].y = cloudlets[num].y;
			}
			bool flag = true;
			while (flag)
			{
				
				AssignmentNEW(V, C, t); 
				
				/* mapping the result of Assignment into R */
				for (int a = 0;a<cloudlet_total;a++)
					for (int b = 0;b<user_total;b++)
						if (S[a][b] == 1)
						{
							R[b][t] = a;

							C[a].b -= users[b][t].b;
							C[a].p -= users[b][t].p;

						}
				flag = false;
				for (int i = 0;i<user_total;i++)
					if (R[i][t] == -1) /* if there is an unassigned user */
					{
						flag = true;
					}

			}

		}

		else
		{
			for (int i = 0;i < user_total;i++)
			{
				R[i][t] = R[i][t - 1];
			}
		}

		/* TEST */
		bool* testarray = new bool[cloudlet_total];
		for (int i = 0;i<cloudlet_total;i++)
			testarray[i] = 0;

		for (int i = 0;i<user_total;i++)
		{
			int index;
			index = R[i][t];
			testarray[index] = 1;

			UsageP[index][t]+=users[i][t].p;
			UsageB[index][t] += users[i][t].b;
		}
		int count = 0;
		for (int i = 0;i<cloudlet_total;i++)
			if (testarray[i] == 1)
				count++;

		load_balance += count;

		/* calculating delay */
		for (int i = 0;i < user_total;i++)
		{
			int index = R[i][t];
			D += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;
			Latency[i] += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;

			/* calculate metrics */
			OffloadingTime += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * theta + (users[i][t].id / users[i][t].b) * 1000000;
			ComputationTime += (users[i][t].w / users[i][t].p) * 1000000;


		}

		if (t >= 2)
			for (int i = 0;i < user_total;i++)
			{

				int index1 = R[i][t - 1];
				int index2 = R[i][t];

				if (index1 != index2)
				{
				    if(i<user_total/3)
					  MigrationNewYork++;
					if(i>=user_total/3 && i<(user_total*2)/3)
					  MigrationOrlando++;
					if(i>= (user_total * 2) / 3 && i<user_total)
					  MigrationDACT++;

					Migrations[i]++;
					migration++;  D += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;
					Latency[i] += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;

					MigrationTime += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;
				}
			}
	}

	return D;
}

void BFD_V2(int** V, int* total, Cloudlets* C, int t)
{

	

	for (int num = 0;num<cloudlet_total;num++)
		C[num] = cloudlets[num];

	for (int j = 0;j < cloudlet_total;j++)
	{
		
			Users* u = new Users[user_total];
			for (int i = 0;i<user_total;i++)
				{
				  users[i][t].index = i;
				  u[i] = users[i][t];
			
				  u[i].total = (u[i].max_p/C[j].p + u[i].max_b/C[j].b)/V[i][j];
				  
				}
				
			std::sort(u, u + user_total, compareTwoApplicationsBFD);

			for (int i = 0;i < user_total;i++)
			{
			  
				if (C[j].p >= u[i].max_p && C[j].b >= u[i].max_b && R[u[i].index][t] ==-1)
				{

					R[u[i].index][t] = j;
					C[j].b -= u[i].max_b;
					C[j].p -= u[i].max_b;
				}
			}
			 
	}

}

void GOAMC(int** V, int* total, Cloudlets* C, int t)
{


	for (int num = 0;num<cloudlet_total;num++)
		{ 
		   C[num] = cloudlets[num];
		   C[num].index = num;
		}

	std::sort(C, C + cloudlet_total, compareTwoCloudletsBFD);

	for (int j = 0;j < cloudlet_total;j++)
	{
		Users* u = new Users[user_total];
		for (int i = 0;i<user_total;i++)
		{
			u[i] = users[i][t];
			u[i].total = (u[i].max_p / C[j].p + u[i].max_b / C[j].b) / V[i][j];

		}

		std::sort(u, u + user_total, compareTwoApplicationsBFD);

		for (int i = 0;i < user_total;i++)
		{
			if (C[j].p >= u[i].max_p && C[j].b >= u[i].max_b && R[u[i].index][t] == -1)
			{
				R[u[i].index][t] = C[j].index;
				C[j].b -= u[i].max_b;
				C[j].p -= u[i].max_b;
			}
		}

	}
}


double Prediction_BFD()
{
	int* total = new int[user_total];
	double D = 0; /* measures delay */

	for (int k = 0;k<user_total;k++)
		if (count_time[k] <= T)
			T = count_time[k];
	T = 20; 
	
	for (int t = 0;t <T;t++)
	{

			for (int i = 0;i < user_total;i++)
			{
				R[i][t] = -1; /* initialization*/

				Users* u = PREDICT_NEW(i, t); /* gets w predicted specifications including u[0], u[1], ..., u[w-1]*/

				total[i] = users[i][t].b + users[i][t].p;

				int max_p = users[i][t].p;
				int max_b = users[i][t].b;
				for (int counter = 0;counter < w;counter++)
				{
					total[i] += (u[counter].b + u[counter].p);
					if (u[counter].b >= max_b)
						max_b = u[counter].b;
					if (u[counter].p >= max_p)
						max_p = u[counter].p;
				}
				users[i][t].total = total[i];
				users[i][t].max_b = max_b;
				users[i][t].max_p = max_p;

				for (int j = 0;j < cloudlet_total;j++) 	/* calculating v_{ij}^t */
				{
					V[i][j] = 0.001*(haversine(users[i][t].x, users[i][t].y, cloudlets[j].x, cloudlets[j].y) * 1000 * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000);
					
					for (int counter = 0;counter < w;counter++)
					{
						V[i][j] += 0.001*(pow(gamma, counter + 1)*(haversine(u[counter].x, u[counter].y, cloudlets[j].x, cloudlets[j].y) * 1000 * theta + (u[counter].id / u[counter].b + u[counter].w / u[counter].p)) * 1000000);
						
					}


					/* checking migration*/
					if (t >= 1 && j != R[i][t - 1])
					{
						int a = R[i][t - 1];
						V[i][j] += 0.001*(haversine(cloudlets[a].x, cloudlets[a].y, cloudlets[j].x, cloudlets[j].y) * 1000 * theta + (users[i][t].md / users[i][t].b) * 1000000);
						
					}

					if (V[i][j]<0)
						V[i][j] = INT_MAX;
				}


			}
			for (int i = 0;i<user_total; i++)
				for (int j = 0;j < cloudlet_total;j++)
					V[i][j] = V[i][j] / 500;

			for (int num = 0;num<cloudlet_total;num++)
			{
				C[num].b = cloudlets[num].b;
				C[num].p = cloudlets[num].p;
				C[num].x = cloudlets[num].x;
				C[num].y = cloudlets[num].y;
			}

			BFD_V2(V, total, C, t);

			
		/* calculating delay */
		for (int i = 0;i < user_total;i++)
		{
			int index = R[i][t];
			D += haversine(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * 1000 * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;
			Latency[i] += haversine(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * 1000 * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;
		}

		if (t >= 1)
			for (int i = 0;i < user_total;i++)
			{

				int index1 = R[i][t - 1];
				int index2 = R[i][t];

				if (index1 != index2)
				{
					Migrations[i]++;
					migration++;  D += haversine(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * 1000 * theta + (users[i][t].md / users[i][t].b) * 1000000;
					Latency[i] += haversine(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * 1000 * theta + (users[i][t].md / users[i][t].b) * 1000000;
				}
			}
	}

	return D;
}

double G-OAMC()
{
	int* total = new int[user_total];
	double D = 0; /* measures delay */


	T = 20; 
	for (int t = 1;t <=T;t++)
	{
		if (t%w == 1 || w==1)
		{


		for (int i = 0;i < user_total;i++)
		{
			R[i][t] = -1; /* initialization*/
						 
			Users* u = PREDICT_NEW(i, t); /* gets w predicted specifications including u[0], u[1], ..., u[w-1]*/
						
			total[i] = users[i][t].b + users[i][t].p;

			int max_p = users[i][t].p;
			int max_b = users[i][t].b;

			int w2; /* w2 is minimum between w and 20-t*/
			if (20 - t<w)
				w2 = 20 - t;
			else
				w2 = w;
			for (int counter = 0;counter < w2;counter++)
			{
				total[i] += (u[counter].b + u[counter].p);
				if (u[counter].b >= max_b)
					max_b = u[counter].b;
				if (u[counter].p >= max_p)
					max_p = u[counter].p;
			}
			users[i][t].total = total[i];
			users[i][t].max_b = max_b;
			users[i][t].max_p = max_p;
			
			for (int j = 0;j < cloudlet_total;j++) 	/* calculating v_{ij}^t */
			{
				V[i][j] = 0.001*(distance(users[i][t].x, users[i][t].y, cloudlets[j].x, cloudlets[j].y) * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000);
				
				int w1; /* w1 is minimum between w and 20-t*/
				if (20 - t<w)
					w1 = 20 - t;
				else
					w1 = w;
				for (int counter = 0;counter < w1;counter++)
				{
					V[i][j] += 0.001*(pow(gamma, counter + 1)*(distance(u[counter].x, u[counter].y, cloudlets[j].x, cloudlets[j].y) * theta + (u[counter].id / u[counter].b + u[counter].w / u[counter].p)) * 1000000);
				
				}

				/* checking migration*/
				if (t >= 2 && j != R[i][t - 1])
				{
					int a = R[i][t - 1];
					V[i][j] += 0.001*(distance(cloudlets[a].x, cloudlets[a].y, cloudlets[j].x, cloudlets[j].y) * theta + (users[i][t].md / users[i][t].b) * 1000000);
				

				}

				if (V[i][j]<0)
					V[i][j] = INT_MAX;
			}


		}
		for (int i = 0;i<user_total; i++)
			for (int j = 0;j < cloudlet_total;j++)
				V[i][j] = V[i][j] / 500;
		for (int num = 0;num<cloudlet_total;num++)
		{
			C[num].b = cloudlets[num].b;
			C[num].p = cloudlets[num].p;
			C[num].x = cloudlets[num].x;
			C[num].y = cloudlets[num].y;
		}
		BFD_V2(V,total, C, t);
		

		}

		else
		{
		 for (int i = 0;i < user_total;i++)
		 {
		 R[i][t] = R[i][t - 1];
		 }
		}
		/* TEST */
		bool* testarray = new bool[cloudlet_total];
		for (int i = 0;i<cloudlet_total;i++)
			testarray[i] = 0;

		for (int i = 0;i<user_total;i++)
		{
			int index;
			index = R[i][t];
			testarray[index] = 1;

			UsageP[index][t] += users[i][t].p;
			UsageB[index][t] += users[i][t].b;
		}
		int count = 0;
		for (int i = 0;i<cloudlet_total;i++)
			if (testarray[i] == 1)
				count++;

		load_balance += count;

		/* calculating delay */
		for (int i = 0;i < user_total;i++)
		{
			int index = R[i][t];

			D += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;
			Latency[i] += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;
			/* calculate metrics */
			OffloadingTime += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * theta + (users[i][t].id / users[i][t].b) * 1000000;
			ComputationTime += (users[i][t].w / users[i][t].p) * 1000000;
		}

		if (t >= 2)
			for (int i = 0;i < user_total;i++)
			{

				int index1 = R[i][t - 1];
				int index2 = R[i][t];

				if (index1 != index2)
				{
					if (i<user_total / 3)
						MigrationNewYork++;
					if (i >= user_total / 3 && i<(user_total * 2) / 3)
						MigrationOrlando++;
					if (i >= (user_total * 2) / 3 && i<user_total)
						MigrationDACT++;

					Migrations[i]++;
					migration++;
					D += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;
					Latency[i] += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;
					MigrationTime += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;
				}
			}
	}


	return D;
}


double Prediction_BFD3()
{
	int* total = new int[user_total];
	double D = 0; /* measures delay */

	for (int k = 0;k<user_total;k++)
		if (count_time[k] <= T)
			T = count_time[k];
	T = 20; 
	for (int t = 0;t <T;t++)
	{
			for (int i = 0;i < user_total;i++)
			{
				R[i][t] = -1; /* initialization*/
							 
				total[i] = users[i][t].b + users[i][t].p;

				int max_p = users[i][t].p;
				int max_b = users[i][t].b;

				users[i][t].total = total[i];
				users[i][t].max_b = max_b;
				users[i][t].max_p = max_p;

				for (int j = 0;j < cloudlet_total;j++) 	/* calculating v_{ij}^t */
				{
					V[i][j] = 0.001*(haversine(users[i][t].x, users[i][t].y, cloudlets[j].x, cloudlets[j].y) * 1000 * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000);


					/* checking migration*/
					if (t >= 1 && j != R[i][t - 1])
					{
						int a = R[i][t - 1];
						V[i][j] += 0.001*(haversine(cloudlets[a].x, cloudlets[a].y, cloudlets[j].x, cloudlets[j].y) * 1000 * theta + (users[i][t].md / users[i][t].b) * 1000000);

					}

					if (V[i][j]<0)
						V[i][j] = INT_MAX;
				}

			}
			for (int i = 0;i<user_total; i++)
				for (int j = 0;j < cloudlet_total;j++)
					V[i][j] = V[i][j] / 500;
		
			for (int num = 0;num<cloudlet_total;num++)
			{
				C[num].b = cloudlets[num].b;
				C[num].p = cloudlets[num].p;
				C[num].x = cloudlets[num].x;
				C[num].y = cloudlets[num].y;
			}

			BFD_V2(V, total, C, t);

		/* calculating delay */
		for (int i = 0;i < user_total;i++)
		{
			int index = R[i][t];
			D += haversine(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * 1000 * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;
			Latency[i] += haversine(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * 1000 * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;
		}

		if (t >= 1)
			for (int i = 0;i < user_total;i++)
			{

				int index1 = R[i][t - 1];
				int index2 = R[i][t];

				if (index1 != index2)
				{
					Migrations[i]++;
					migration++;  D += haversine(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * 1000 * theta + (users[i][t].md / users[i][t].b) * 1000000;
					Latency[i] += haversine(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * 1000 * theta + (users[i][t].md / users[i][t].b) * 1000000;
				}
			}
	}

	return D;
}



void IP()
{
	T = 20;
	IloEnv env;
	try {
		IloModel model(env);


		/* define the bool vars here for the 4-D matrix */
		BoolVar4Matrix alpha(env, user_total);
		/* initialize this matrix */
		for (int i = 0; i< user_total; i++) {
			alpha[i] = BoolVar3Matrix(env, cloudlet_total);
			for (int j = 0; j< cloudlet_total; j++) {
				alpha[i][j] = BoolVarMatrix(env, cloudlet_total);
				for (int k = 0; k<cloudlet_total; k++) {
					alpha[i][j][k] = IloBoolVarArray(env, T+1);
					for (int t = 1; t<=T; t++) {
						alpha[i][j][k][t] = IloBoolVar(env);

					}
				}

			}
		}

		/* define the bool vars here for the 3-D matrix */
		BoolVar3Matrix mu(env, user_total);
		/* initialize this matrix */
		for (int i = 0; i< user_total; i++) {
			mu[i] = BoolVarMatrix(env, cloudlet_total);
			for (int j = 0; j< cloudlet_total; j++) {
				mu[i][j] = IloBoolVarArray(env, T+1);
				for (int t = 1; t<=T; t++) {
					mu[i][j][t] = IloBoolVar(env);
				}
			}
		}

		IloExpr expression1(env);
		IloExpr expression2(env);

		for (int i = 0;i < user_total;i++)
		{
			for (int j = 0;j < cloudlet_total;j++)
			{
				for (int t = 1;t <=T;t++)
				{
					expression1 += (distance(users[i][t].x, users[i][t].y, cloudlets[j].x, cloudlets[j].y) * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000)*mu[i][j][t];
				}
			}
		}

		for (int i = 0;i < user_total;i++)
		{
			for (int j = 0;j < cloudlet_total;j++)
			{
				for (int k = 0;k < cloudlet_total;k++)
				{
					for (int t = 2;t <=T;t++)
					{
						if (j != k)
							expression2 += (distance(cloudlets[j].x, cloudlets[j].y, cloudlets[k].x, cloudlets[k].y) * theta + (users[i][t].md / users[i][t].b) * 1000000)*alpha[i][j][k][t];
					}
				}

			}
		}

		model.add(IloMinimize(env, expression1 + expression2)); /* minimize the cost */

																/* constraints */

		IloExpr Expr(env);

		for (int i = 0;i < user_total;i++) /* 1 */
		{
			for (int t = 1;t <=T;t++)
			{
				for (int j = 0;j < cloudlet_total;j++)
				{
					Expr += mu[i][j][t];
				}

				model.add(Expr == 1);
				Expr.end();
				Expr = IloExpr(env);
			}
		}

		IloExpr Expr2(env);

		for (int j = 0;j < cloudlet_total;j++)
		{
			for (int t = 1;t <=T;t++)
			{
				for (int i = 0;i < user_total;i++)
				{
					Expr2 += users[i][t].p * mu[i][j][t];
				}

				model.add(Expr2 <= cloudlets[j].p);
				Expr2.end();
				Expr2 = IloExpr(env);
			}
		}

		IloExpr Expr3(env);

		for (int j = 0;j < cloudlet_total;j++)
		{
			for (int t = 1;t <=T;t++)
			{
				for (int i = 0;i < user_total;i++)
				{
					Expr3 += users[i][t].b * mu[i][j][t];
				}

				model.add(Expr3 <= cloudlets[j].b);
				Expr3.end();
				Expr3 = IloExpr(env);
			}
		}

		for (int i = 0;i < user_total;i++)
		{
			for (int j = 0;j < cloudlet_total;j++)
			{
				for (int k = 0;k < cloudlet_total;k++)
				{
					for (int t = 2;t <=T;t++)
					{
						if (j != k)
							model.add(alpha[i][j][k][t] >= mu[i][j][t - 1] + mu[i][k][t] - 1);
					}
				}
			}
		}
		IloCplex cplex(env);
		//cplex.setParam(IloCplex::EpGap,0.02);
		cplex.setParam(IloCplex::Param::TimeLimit, 3600);

		cplex.extract(model);
		cplex.exportModel("model.lp");
		env.out() << "Variables binarias: " << cplex.getNbinVars() << endl;
		env.out() << "Variables Enteras: " << cplex.getNintVars() << endl;
		env.out() << "Filas - Restricciones: " << cplex.getNrows() << endl;
		env.out() << "Columnas - Variables: " << cplex.getNcols() << endl;

		if (!cplex.solve()) {
			env.error() << "No es posible resolver :-(" << endl;
			throw (-1);
		}

		env.out() << "Es optimo ? = " << cplex.getStatus() << endl;
		env.out() << "Valor de fo = " << cplex.getObjValue() << endl;
		env.out() << "Se demoro   = " << env.getTime() << endl;

		
		/* this for loop is just for calculating load balancing and P utilization */
		for (int t = 1;t <= T;t++)
		{
			/* TEST */
			bool* testarray = new bool[cloudlet_total];
			for (int i = 0;i<cloudlet_total;i++)
				testarray[i] = 0;

			for (int i = 0;i<user_total;i++)
				for (int j = 0;j < cloudlet_total;j++)
					if (cplex.getValue(mu[i][j][t]) == 1)
					{
						testarray[j] = 1;

						UsageP[j][t] += users[i][t].p;
						UsageB[j][t] += users[i][t].b;
					}
			int count = 0;
			for (int i = 0;i<cloudlet_total;i++)
				if (testarray[i] == 1)
					count++;

			load_balance += count;
		}
		 
		

		int migration2 = 0;
		for (int i = 0;i < user_total;i++)
		{
			for (int j = 0;j < cloudlet_total;j++)
			{
				for (int k = 0;k < cloudlet_total;k++)
				{
					for (int t = 2;t <=T;t++)
					{

						if (j != k)
							if (cplex.getValue(alpha[i][j][k][t]) == 1)
							{
								if (i<user_total / 3)
									MigrationNewYork++;
								if (i >= user_total / 3 && i<(user_total * 2) / 3)
									MigrationOrlando++;
								if (i >= (user_total * 2) / 3 && i<user_total)
									MigrationDACT++;

								 migration2++;
								 Migrations[i]++; 					
								 Latency[i] += distance(cloudlets[j].x, cloudlets[j].y, cloudlets[k].x, cloudlets[k].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;
								 MigrationTime += distance(cloudlets[j].x, cloudlets[j].y, cloudlets[k].x, cloudlets[k].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;

							}
					}
				}
			}
		}

		for (int i = 0;i < user_total;i++)
		{
			for (int j = 0;j < cloudlet_total;j++)
			{
				
					for (int t = 1;t <=T;t++)
					{

							if (cplex.getValue(mu[i][j][t]) == 1)
							{
							 Latency[i]+= distance(users[i][t].x, users[i][t].y, cloudlets[j].x, cloudlets[j].y) * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;

							 /* calculate metrics */
							 OffloadingTime += distance(users[i][t].x, users[i][t].y, cloudlets[j].x, cloudlets[j].y) * theta + (users[i][t].id / users[i][t].b) * 1000000;
							 ComputationTime += (users[i][t].w / users[i][t].p) * 1000000;
							}
					}
			
			}
		}

		cout << "IP Migration " << migration2 << endl;
	}
	catch (IloException& ex) {
		cerr << "Error Cplex: " << ex << endl;
	}
	catch (...) {
		cerr << "Error Cpp" << endl;
	}
	env.end();

}

bool compareTwoApplicationsForBFD(Users user_1, Users user_2)
{
	return user_1.p>=user_2.p &&  user_1.b >= user_2.b;
}


double BFD()
{

	
	double D = 0;


	T = 20; 
	for (int t = 1;t <=T;t++)
	{
		Users* u = new Users[user_total];
		for (int i = 0;i<user_total;i++)
		{
		     users[i][t].index = i;
			 u[i] = users[i][t];
		}
		std::sort(u, u + user_total, compareTwoApplicationsForBFD);
	

		for (int num = 0;num<cloudlet_total;num++)
			C[num] = cloudlets[num];

		for (int i = 0;i < user_total;i++)
		{
			double MAX = INT_MIN;
			int index = -1;
			for (int j = 0;j < cloudlet_total;j++)
			{
			    
				double fp = (u[i].p ) / (C[j].p + 1);
				double fb = (u[i].b) / (C[j].b + 1);

				if (fp + fb >= MAX && C[j].p >= u[i].p && C[j].b >= u[i].b)
				{
					index = j;
					MAX = fp + fb;
				}
			}
			

			R[u[i].index][t] = index;
			C[index].b -= u[i].b;
			C[index].p -= u[i].p;
		}

		/* TEST */
		bool* testarray = new bool[cloudlet_total];
		for (int i = 0;i<cloudlet_total;i++)
			testarray[i] = 0;

		for (int i = 0;i<user_total;i++)
		{
			int index;
			index = R[i][t];
			testarray[index] = 1;

			UsageP[index][t] += users[i][t].p;
			UsageB[index][t] += users[i][t].b;

		}
		int count = 0;
		for (int i = 0;i<cloudlet_total;i++)
			if (testarray[i] == 1)
				count++;

		load_balance += count;

		for (int i = 0;i < user_total;i++)
		{
			int index = R[i][t];
			D += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y)  * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;
			Latency[i] += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;

			/* calculate metrics */
			OffloadingTime += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * theta + (users[i][t].id / users[i][t].b) * 1000000;
			ComputationTime += (users[i][t].w / users[i][t].p) * 1000000;
		}
		if (t >= 2)
			for (int i = 0;i < user_total;i++)
			{

				int index1 = R[i][t - 1];
				int index2 = R[i][t];

				if (index1 != index2)
				{

					if (i<user_total / 3)
						MigrationNewYork++;
					if (i >= user_total / 3 && i<(user_total * 2) / 3)
						MigrationOrlando++;
					if (i >= (user_total * 2) / 3 && i<user_total)
						MigrationDACT++;

				    Migrations[i]++;
					migration++; D += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;
					Latency[i] += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;
					MigrationTime += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;

				}
			}
	}
	return D;
}
bool* temporary1()
{
	bool s[100];
	for (int i = 0;i<100;i++)
		s[i] = 1;

	Struct a;
	for (int i = 0;i<100;i++)
		a.temp[i] = 1;
	return a.temp;
}

double average = 0;
int num = 0;
void JustTest(int t)
{
	/* for predicting x */
	for (int i = 0;i < user_total;i++)
	{
		for (int h = 0;h < w;h++)
		{
			/* for predicting x */
			temporary[h].x = users[i][t].x;
			double dif = users[i][t].x - users[i][1].x;

			if (t != 0)
				temporary[h].x += (double(h + 1) / (double)t)*dif;

			/* for predicting y */
			temporary[h].y = users[i][t].y;
			double dif1 = users[i][t].y - users[i][1].y;

			if (t != 0)
				temporary[h].y += (double(h + 1) / (double)t)*dif1;
	   }
	   cout << "users[i][t].x " << users[i][t].x <<" users[i][t].y "<< users[i][t].y<<endl;

	   cout<<"  users[i][t+1].x "<< users[i][t+1].x<<" predicted "<<  temporary[0].x<< "percentage of closeness "<< abs(users[i][t + 1].x- temporary[0].x)/(max_x-min_x)<<endl;
	   cout<< " users[i][t+2].x " << users[i][t + 2].x << " predicted " << temporary[1].x << "percentage of closeness "<< abs(users[i][t + 2].x - temporary[1].x) /(max_x - min_x) << endl;
	   cout<< " users[i][t+3].x " << users[i][t + 3].x << " predicted " << temporary[2].x << "percentage of closeness "<< abs(users[i][t + 3].x - temporary[2].x) /(max_x - min_x) << endl;
	   cout<< " users[i][t+4].x " << users[i][t + 4].x << " predicted " << temporary[3].x << "percentage of closeness "<< abs(users[i][t + 4].x - temporary[3].x) /(max_x - min_x) << endl;
	   cout<< " users[i][t+5].x " << users[i][t + 5].x << " predicted " << temporary[4].x << "percentage of closeness "<<abs(users[i][t + 5].x - temporary[4].x) /(max_x - min_x) << endl;

	   cout << " users[i][t+1].y " << users[i][t + 1].y << " predicted " << temporary[0].y << "percentage of closeness "<<abs(users[i][t + 1].y - temporary[0].y) /(max_y - min_y) << endl;
	   cout << " users[i][t+2].y " << users[i][t + 2].y << " predicted " << temporary[1].y <<"percentage of closeness "<< abs(users[i][t + 2].y - temporary[1].y) /(max_y - min_y) << endl;
	   cout << " users[i][t+3].y " << users[i][t + 3].y << " predicted " << temporary[2].y <<"percentage of closeness "<< abs(users[i][t + 3].y - temporary[2].y) /(max_y - min_y) << endl;
	   cout << " users[i][t+4].y " << users[i][t + 4].y << " predicted " << temporary[3].y <<"percentage of closeness "<< abs(users[i][t + 4].y - temporary[3].y) /(max_y - min_y) << endl;
	   cout << " users[i][t+5].y " << users[i][t + 5].y << " predicted " << temporary[4].y <<"percentage of closeness "<< abs(users[i][t + 5].y - temporary[4].y) /(max_y - min_y) << endl;
	   cout<<"for t+1 "<< haversine(temporary[0].x, temporary[0].y, users[i][t+1].x, users[i][t+1].y) /haversine(users[i][t + 1].x, users[i][t + 1].y, users[i][t].x, users[i][t].y) << endl;

	   
	   cout<<"first "<< haversine(users[i][t].x, users[i][t].y, users[i][t + 1].x, users[i][t + 1].y)<<endl;

	   cout << "second " << haversine(users[i][t+1].x, users[i][t+1].y, users[i][t + 2].x, users[i][t + 2].y) << endl;

	   cout << "third " << haversine(users[i][t+2].x, users[i][t+2].y, users[i][t + 3].x, users[i][t + 3].y) << endl;

	   cout << "fourth " << haversine(users[i][t+3].x, users[i][t+3].y, users[i][t + 4].x, users[i][t + 4].y) << endl;

	   cout << "fifth " << haversine(users[i][t+4].x, users[i][t+4].y, users[i][t + 5].x, users[i][t + 5].y) << endl;

	   average += haversine(temporary[0].x, temporary[0].y, users[i][t + 1].x, users[i][t + 1].y) / (0.0001+haversine(users[i][t].x, users[i][t].y,users[i][t + 1].x, users[i][t + 1].y));
	   num++;
	 }

}

void NewDataSet()
{
	/* new dataset --- beginning-----  */
	for (int i = 0;i<user_total;i++)
	{
		users[i] = new Users[200];
		R[i] = new int[200];
		count_time[i] = 200;
	}


	ifstream data("C:\\Users\\Erfan\\Desktop\\DACT Strict-Dataset.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';

	string first;
	string second;
	string third;
	string fourth;
	string fifth;
	string sixth;
	string seventh;
	string eighth;
	string ninth;
	string tenth;
	string eleventh;
	string word;
	while (data.good())
	{

		getline(data, first, ',');
		getline(data, second, ',');
		getline(data, third, ',');
		getline(data, fourth, ',');
		getline(data, fifth, ',');
		getline(data, sixth, ',');
		getline(data, seventh, ',');
		getline(data, eighth, ',');
		getline(data, ninth, ',');
		getline(data, tenth, ',');
		

		stringstream geek(eighth);
		double x = 0;
		geek >> x;
		cout.precision(eighth.size());


		if (x <= min_x)
			min_x = x;
		if (x>max_x)
			max_x = x;

		stringstream geek2(ninth);
		double y = 0;
		geek2 >> y;
		cout.precision(ninth.size());

		if (y <= min_y)
			min_y = y;
		if (y>max_y)
			max_y = y;

	}
	/* end */
	ifstream data1("C:\\Users\\Erfan\\Desktop\\DACT Strict-Dataset.csv");
     
	int count = 0;
	int user_count = 0;
	int index = 0;

	if (!data1.is_open()) cout << "ERROR: File Open" << '\n';
	while (data1.good() && user_count<user_total)
	{

		getline(data1, first, ',');
		getline(data1, second, ',');
		getline(data1, third, ',');
		getline(data1, fourth, ',');
		getline(data1, fifth, ',');
		getline(data1, sixth, ',');
		getline(data1, seventh, ',');
		getline(data1, eighth, ',');
		getline(data1, ninth, ',');
		getline(data1, tenth, ',');
	
		count++;
		stringstream geek(eighth);
		double x = 0;
		geek >> x;
		cout.precision(eighth.size());


		if (x <= min_x)
			min_x = x;
		if (x>max_x)
			max_x = x;

		stringstream geek2(ninth);
		double y = 0;
		geek2 >> y;
		cout.precision(ninth.size());


		if (y <= min_y)
			min_y = y;
		if (y>max_y)
			max_y = y;

		if(count%10==1)    
		{
			users[user_count][index].x = x;
			users[user_count][index].y = y;
			index++;
		}
		if (index == 21)
		{
			index = 0;
			user_count++;
			count = 0;
		}

	}

}

void SpecPredict()
{
	for (int i = 0;i<user_total;i++)
	{
		temporary_test[i] = new Users[20];
	}

	for (int t = 1;t <= 19;t++)
	{
		for (int i = 0;i < user_total;i++)
		{
			for (int h = 0;h < 20 - t;h++)
			{
				temporary_test[i][h].x = users[i][t].x;
				double dif = users[i][t].x - users[i][1].x;
				if (t != 0)
					temporary_test[i][h].x += (double(h + 1) / double(t))*dif;
				temporary_test[i][h].x = (int)temporary_test[i][h].x;

				/* for predicting id */
				double sum1 = 0;
				for (int counter = 1; counter <= t;counter++)
					sum1 += (counter)*users[i][counter].id;
				sum1 /= (t)*(t + 1) / 2;
				temporary_test[i][h].id = sum1;

				/* for predicting y */
				temporary_test[i][h].y = users[i][t].y;
				double dif1 = users[i][t].y - users[i][1].y;
				if (t != 0)
					temporary_test[i][h].y += (double(h + 1) / double(t))*dif1;
				temporary_test[i][h].y = (int)temporary_test[i][h].y;
				
				/* for predicting md */
				int sum2 = 0;
				for (int counter = 1; counter <= t;counter++)
					sum2 += (counter)*users[i][counter].md;
				sum2 /= ((t)*(t + 1) / 2);
				temporary_test[i][h].md = sum2;

				/* for predicting w */
				int sum3 = 0;
				for (int counter = 1; counter <= t;counter++)
					sum3 += (counter)*users[i][counter].w;
				sum3 /= ((t)*(t + 1) / 2);
				temporary_test[i][h].w = sum3;

				/* for predicting p */
				int sum4 = 0;
				for (int counter = 1; counter <= t;counter++)
					sum4 += (counter)*users[i][counter].p;
				sum4 /= ((t)*(t + 1) / 2);
				temporary_test[i][h].p = sum4;

				/* for predicting b */
				int sum5 = 0;
				for (int counter = 1; counter <= t;counter++)
					sum5 += (counter)*users[i][counter].b;
				sum5 /= ((t)*(t + 1) / 2);
				temporary_test[i][h].b = sum5;

		    }
		}

		/* write to files */
		/*------------------------------------------*/
		stringstream out;
		out << t;
		ofstream myfile;
		myfile.open("C:\\Users\\Erfan\\Desktop\\Results\\x\\Results_50\\results_Sample_Data_50_x_"+ out.str() + ".csv");
		myfile.precision(15);
		for (int i = 0;i < user_total;i++)
		{
		   for(int n=1;n<=t;n++)
		      myfile << users[i][n].x<<",";
		   for(int h=0;h<20-t;h++)
			myfile << temporary_test[i][h].x << "," ;

		   myfile<<endl;
		}

		/*---------------------------------------*/
		stringstream out2;
		out2 << t;
		ofstream myfile2;
		myfile2.open("C:\\Users\\Erfan\\Desktop\\Results\\y\\Results_50\\results_Sample_Data_50_y_" + out2.str() + ".csv");
		myfile2.precision(15);
		for (int i = 0;i < user_total;i++)
		{
			for (int n = 1;n <= t;n++)
				myfile2 << users[i][n].y << ",";
			for (int h = 0;h<20 - t;h++)
				myfile2 << temporary_test[i][h].y << ",";

			myfile2 << endl;
		}

		/*---------------------------------------*/
		stringstream out3;
		out3 << t;
		ofstream myfile3;
		myfile3.open("C:\\Users\\Erfan\\Desktop\\Results\\md\\Results_50\\results_Sample_Data_50_md_" + out3.str() + ".csv");
		myfile3.precision(15);
		for (int i = 0;i < user_total;i++)
		{
			for (int n = 1;n <= t;n++)
				myfile3 << users[i][n].md << ",";
			for (int h = 0;h<20 - t;h++)
				myfile3 << temporary_test[i][h].md << ",";

			myfile3 << endl;
		}

		/*---------------------------------------*/
		stringstream out4;
		out4 << t;
		ofstream myfile4;
		myfile4.open("C:\\Users\\Erfan\\Desktop\\Results\\id\\Results_50\\results_Sample_Data_50_id_" + out4.str() + ".csv");
		myfile4.precision(15);
		for (int i = 0;i < user_total;i++)
		{
			for (int n = 1;n <= t;n++)
				myfile4 << users[i][n].id << ",";
			for (int h = 0;h<20 - t;h++)
				myfile4 << temporary_test[i][h].id << ",";

			myfile4 << endl;
		}

		/*---------------------------------------*/
		stringstream out5;
		out5 << t;
		ofstream myfile5;
		myfile5.open("C:\\Users\\Erfan\\Desktop\\Results\\w\\Results_50\\results_Sample_Data_50_w_" + out5.str() + ".csv");
		myfile5.precision(15);
		for (int i = 0;i < user_total;i++)
		{
			for (int n = 1;n <= t;n++)
				myfile5 << users[i][n].w << ",";
			for (int h = 0;h<20 - t;h++)
				myfile5 << temporary_test[i][h].w << ",";

			myfile5 << endl;
		}

		/*---------------------------------------*/
		stringstream out6;
		out6 << t;
		ofstream myfile6;
		myfile6.open("C:\\Users\\Erfan\\Desktop\\Results\\p\\Results_50\\results_Sample_Data_50_p_" + out6.str() + ".csv");
		myfile6.precision(15);
		for (int i = 0;i < user_total;i++)
		{
			for (int n = 1;n <= t;n++)
				myfile6 << users[i][n].p << ",";
			for (int h = 0;h<20 - t;h++)
				myfile6 << temporary_test[i][h].p << ",";

			myfile6 << endl;
		}

		/*---------------------------------------*/
		stringstream out7;
		out7 << t;
		ofstream myfile7;
		myfile7.open("C:\\Users\\Erfan\\Desktop\\Results\\b\\Results_50\\results_Sample_Data_50_b_" + out7.str() + ".csv");
		myfile7.precision(15);
		for (int i = 0;i < user_total;i++)
		{
			for (int n = 1;n <= t;n++)
				myfile7 << users[i][n].b << ",";
			for (int h = 0;h<20 - t;h++)
				myfile7 << temporary_test[i][h].b << ",";

			myfile7 << endl;
		}

	}

}
void Read()
{

	/* new dataset --- beginning-----  */
	for (int i = 0;i<user_total;i++)
	{
		users[i] = new Users[200];
		R[i] = new int[200];
		count_time[i] = 200;
	}


	cout.precision(15);

	ifstream data("C:\\Users\\Erfan\\Desktop\\Datasets\\Sample Data_50.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';

	string first;
	string second;
	string third;
	string fourth;
	string fifth;
	string sixth;
	string seventh;
	int count = 0;
	int user_count = 0;
	int index = 1;
	while (data.good() && user_count<user_total)
	{

	    count++;
		getline(data, first, ',');
		getline(data, second, ',');
		getline(data, third, ',');
		getline(data, fourth, ',');
		getline(data, fifth, ',');
		getline(data, sixth, ',');
		getline(data, seventh, '\n');
		
		stringstream geek(first);
		double x = 0;
		geek >> x;

		users[user_count][index].x = x;

		stringstream geek2(second);
		double y = 0;
		geek2 >> y;

		users[user_count][index].y = y;

		const double halfC = M_PI / 180;

		double lat = x;
		if(index == 1 )
		 cout<<lat<<"   ";
		double lon = y;
		double h = 274.93;

		double LatRad = lat*halfC;
		double LongRad = lon*halfC;

		double a2 = 6378.1370 * 1000;
		double b2 = 6356.7523 * 1000;
		double N;
		double e = 1 - pow(b2, 2) / pow(a2, 2);
		N = a2 / sqrt(1.0 - (e*pow(sin(LatRad), 2)));
		double cosLatRad = cos(LatRad);
		double cosLongiRad = cos(LongRad);
		double sinLatRad = sin(LatRad);
		double sinLongiRad = sin(LongRad);

		int x1 = (N + h)*cosLatRad*cosLongiRad;
		int y1 = (N + h)*cosLatRad*sinLongiRad;
		double z = ((pow(b2, 2) / pow(a2, 2))*N + h)*sinLatRad;

		users[user_count][index].x = x1;
		users[user_count][index].y = y1;

		if (index == 1)
			cout << x1 << endl;

		stringstream geek3(third);
		double md = 0;
		geek3 >> md;

		users[user_count][index].md = md;

		stringstream geek4(fourth);
		double id = 0;
		geek4 >> id;

		users[user_count][index].id = id;

		stringstream geek5(fifth);
		double w = 0;
		geek5 >> w;

		users[user_count][index].w = w;

		stringstream geek6(sixth);
		double p = 0;
		geek6 >> p;

		users[user_count][index].p = p;

		stringstream geek7(seventh);
		double b = 0;
		geek7 >> b;

		users[user_count][index].b = b;

		index++;

		if (index == 21)
		{
			getline(data, first, ',');
			getline(data, second, ',');
			getline(data, third, ',');
			getline(data, fourth, ',');
			getline(data, fifth, ',');
			getline(data, sixth, ',');

			index = 1;
			user_count++;
			count = 0;
		}
		
	}
	
	stringstream out6;
	ofstream myfile6;
	myfile6.open("C:\\Users\\Erfan\\Desktop\\Datasets\\Sample_Data_Conversion_50_x.csv");
	myfile6.precision(15);
	for (int i = 0;i < user_total;i++)
	{
		for (int n = 1;n <= 20;n++)
			myfile6 << users[i][n].x << ",";
		myfile6 << endl;
	}

	stringstream out7;
	ofstream myfile7;
	myfile7.open("C:\\Users\\Erfan\\Desktop\\Datasets\\Sample_Data_Conversion_50_y.csv");
	myfile7.precision(15);
	for (int i = 0;i < user_total;i++)
	{
		for (int n = 1;n <= 20;n++)
			myfile7 << users[i][n].y << ",";
		myfile7 << endl;
	}

}

void Results()
{
	/* new dataset --- beginning-----  */
	for (int i = 0;i<user_total;i++)
	{
		users[i] = new Users[200];
		R[i] = new int[200];

		first_prediction[i] = new Users[20];
		second_prediction[i] = new Users[20];
		third_prediction[i] = new Users[20];
		fourth_prediction[i] = new Users[20];
		fifth_prediction[i] = new Users[20];

		count_time[i] = 200;
	}


	cout.precision(15);

	int t = 1;
	while(t<=19)
	{ 
	 stringstream out;
	 out << t;
	
	ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\x\\Results_50\\results_Sample_Data_50_x_" + out.str() + ".csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';

	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;

	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty, '\n');



		/*  ---------------------------   */
		stringstream geek(one);
		double x1 = 0;
		geek >> x1;

		stringstream geek2(two);
		double x2 = 0;
		geek2 >> x2;

		stringstream geek3(three);
		double x3 = 0;
		geek3 >> x3;

		stringstream geek4(four);
		double x4 = 0;
		geek4 >> x4;

		stringstream geek5(five);
		double x5 = 0;
		geek5 >> x5;

		stringstream geek6(six);
		double x6 = 0;
		geek6 >> x6;

		stringstream geek7(seven);
		double x7 = 0;
		geek7 >> x7;

		stringstream geek8(eight);
		double x8 = 0;
		geek8 >> x8;

		stringstream geek9(nine);
		double x9 = 0;
		geek9 >> x9;

		stringstream geek10(ten);
		double x10 = 0;
		geek10 >> x10;

		stringstream geek11(eleven);
		double x11 = 0;
		geek11 >> x11;

		stringstream geek12(twelve);
		double x12 = 0;
		geek12 >> x12;

		stringstream geek13(thirteen);
		double x13 = 0;
		geek13 >> x13;

		stringstream geek14(fourteen);
		double x14 = 0;
		geek14 >> x14;

		stringstream geek15(fifteen);
		double x15 = 0;
		geek15 >> x15;

		stringstream geek16(sixteen);
		double x16 = 0;
		geek16 >> x16;

		stringstream geek17(seventeen);
		double x17 = 0;
		geek17 >> x17;

		stringstream geek18(eighteen);
		double x18 = 0;
		geek18 >> x18;

		stringstream geek19(nineteen);
		double x19 = 0;
		geek19 >> x19;

		stringstream geek20(twenty);
		double x20 = 0;
		geek20 >> x20;

		if (t == 1)
		{
			first_prediction[user_count][t].x = x2;
			second_prediction[user_count][t].x = x3;
			third_prediction[user_count][t].x = x4;
			fourth_prediction[user_count][t].x = x5;
			fifth_prediction[user_count][t].x = x6;

		}
		if (t == 2)
		{
			first_prediction[user_count][t].x = x3;
			second_prediction[user_count][t].x = x4;
			third_prediction[user_count][t].x = x5;
			fourth_prediction[user_count][t].x = x6;
			fifth_prediction[user_count][t].x = x7;

		}
		if (t == 3)
		{
			first_prediction[user_count][t].x = x4;
			second_prediction[user_count][t].x = x5;
			third_prediction[user_count][t].x = x6;
			fourth_prediction[user_count][t].x = x7;
			fifth_prediction[user_count][t].x = x8;

		}
		if (t == 4)
		{
			first_prediction[user_count][t].x = x5;
			second_prediction[user_count][t].x = x6;
			third_prediction[user_count][t].x = x7;
			fourth_prediction[user_count][t].x = x8;
			fifth_prediction[user_count][t].x = x9;

		}
		if (t == 5)
		{
			first_prediction[user_count][t].x = x6;
			second_prediction[user_count][t].x = x7;
			third_prediction[user_count][t].x = x8;
			fourth_prediction[user_count][t].x = x9;
			fifth_prediction[user_count][t].x = x10;

		}
		if (t == 6)
		{
			first_prediction[user_count][t].x = x7;
			second_prediction[user_count][t].x = x8;
			third_prediction[user_count][t].x = x9;
			fourth_prediction[user_count][t].x = x10;
			fifth_prediction[user_count][t].x = x11;

		}
		if (t == 7)
		{
			first_prediction[user_count][t].x = x8;
			second_prediction[user_count][t].x = x9;
			third_prediction[user_count][t].x = x10;
			fourth_prediction[user_count][t].x = x11;
			fifth_prediction[user_count][t].x = x12;

		}
		if (t == 8)
		{
			first_prediction[user_count][t].x = x9;
			second_prediction[user_count][t].x = x10;
			third_prediction[user_count][t].x = x11;
			fourth_prediction[user_count][t].x = x12;
			fifth_prediction[user_count][t].x = x13;

		}
		if (t == 9)
		{
			first_prediction[user_count][t].x = x10;
			second_prediction[user_count][t].x = x11;
			third_prediction[user_count][t].x = x12;
			fourth_prediction[user_count][t].x = x13;
			fifth_prediction[user_count][t].x = x14;

		}
		if (t == 10)
		{
			first_prediction[user_count][t].x = x11;
			second_prediction[user_count][t].x = x12;
			third_prediction[user_count][t].x = x13;
			fourth_prediction[user_count][t].x = x14;
			fifth_prediction[user_count][t].x = x15;

		}
		if (t == 11)
		{
			first_prediction[user_count][t].x = x12;
			second_prediction[user_count][t].x = x13;
			third_prediction[user_count][t].x = x14;
			fourth_prediction[user_count][t].x = x15;
			fifth_prediction[user_count][t].x = x16;

		}
		if (t == 12)
		{
			first_prediction[user_count][t].x = x13;
			second_prediction[user_count][t].x = x14;
			third_prediction[user_count][t].x = x15;
			fourth_prediction[user_count][t].x = x16;
			fifth_prediction[user_count][t].x = x17;

		}
		if (t == 13)
		{
			first_prediction[user_count][t].x = x14;
			second_prediction[user_count][t].x = x15;
			third_prediction[user_count][t].x = x16;
			fourth_prediction[user_count][t].x = x17;
			fifth_prediction[user_count][t].x = x18;

		}
		if (t == 14)
		{
			first_prediction[user_count][t].x = x15;
			second_prediction[user_count][t].x = x16;
			third_prediction[user_count][t].x = x17;
			fourth_prediction[user_count][t].x = x18;
			fifth_prediction[user_count][t].x = x19;

		}
		if (t == 15)
		{
			first_prediction[user_count][t].x = x16;
			second_prediction[user_count][t].x = x17;
			third_prediction[user_count][t].x = x18;
			fourth_prediction[user_count][t].x = x19;
			fifth_prediction[user_count][t].x = x20;

		}
		if (t == 16)
		{
			first_prediction[user_count][t].x = x17;
			second_prediction[user_count][t].x = x18;
			third_prediction[user_count][t].x = x19;
			fourth_prediction[user_count][t].x = x20;

		}
		if (t == 17)
		{
			first_prediction[user_count][t].x = x18;
			second_prediction[user_count][t].x = x19;
			third_prediction[user_count][t].x = x20;

		}
		if (t == 18)
		{
			first_prediction[user_count][t].x = x19;
			second_prediction[user_count][t].x = x20;

		}
		if (t == 19)
		{
			first_prediction[user_count][t].x = x20;

		}


		user_count++;
	 }
	t++;
	}

	/* y */
	{
	int t = 1;
	while (t <= 19)
	{
		stringstream out;
		out << t;
		
		ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\y\\Results_50\\results_Sample_Data_50_y_" + out.str() + ".csv");

		if (!data.is_open()) cout << "ERROR: File Open" << '\n';

		string one;
		string two;
		string three;
		string four;
		string five;
		string six;
		string seven;
		string eight;
		string nine;
		string ten;
		string eleven;
		string twelve;
		string thirteen;
		string fourteen;
		string fifteen;
		string sixteen;
		string seventeen;
		string eighteen;
		string nineteen;
		string twenty;

		int user_count = 0;
		
		while (data.good() && user_count<user_total)
		{

			getline(data, one, ',');
			getline(data, two, ',');
			getline(data, three, ',');
			getline(data, four, ',');
			getline(data, five, ',');
			getline(data, six, ',');
			getline(data, seven, ',');
			getline(data, eight, ',');
			getline(data, nine, ',');
			getline(data, ten, ',');
			getline(data, eleven, ',');
			getline(data, twelve, ',');
			getline(data, thirteen, ',');
			getline(data, fourteen, ',');
			getline(data, fifteen, ',');
			getline(data, sixteen, ',');
			getline(data, seventeen, ',');
			getline(data, eighteen, ',');
			getline(data, nineteen, ',');
			getline(data, twenty, '\n');



			/*  ---------------------------   */
			stringstream geek(one);
			double y1 = 0;
			geek >> y1;

			stringstream geek2(two);
			double y2 = 0;
			geek2 >> y2;

			stringstream geek3(three);
			double y3 = 0;
			geek3 >> y3;

			stringstream geek4(four);
			double y4 = 0;
			geek4 >> y4;

			stringstream geek5(five);
			double y5 = 0;
			geek5 >> y5;

			stringstream geek6(six);
			double y6 = 0;
			geek6 >> y6;

			stringstream geek7(seven);
			double y7 = 0;
			geek7 >> y7;

			stringstream geek8(eight);
			double y8 = 0;
			geek8 >> y8;

			stringstream geek9(nine);
			double y9 = 0;
			geek9 >> y9;

			stringstream geek10(ten);
			double y10 = 0;
			geek10 >> y10;

			stringstream geek11(eleven);
			double y11 = 0;
			geek11 >> y11;

			stringstream geek12(twelve);
			double y12 = 0;
			geek12 >> y12;

			stringstream geek13(thirteen);
			double y13 = 0;
			geek13 >> y13;

			stringstream geek14(fourteen);
			double y14 = 0;
			geek14 >> y14;

			stringstream geek15(fifteen);
			double y15 = 0;
			geek15 >> y15;

			stringstream geek16(sixteen);
			double y16 = 0;
			geek16 >> y16;

			stringstream geek17(seventeen);
			double y17 = 0;
			geek17 >> y17;

			stringstream geek18(eighteen);
			double y18 = 0;
			geek18 >> y18;

			stringstream geek19(nineteen);
			double y19 = 0;
			geek19 >> y19;

			stringstream geek20(twenty);
			double y20 = 0;
			geek20 >> y20;

			if (t == 1)
			{
				first_prediction[user_count][t].y = y2;
				second_prediction[user_count][t].y = y3;
				third_prediction[user_count][t].y = y4;
				fourth_prediction[user_count][t].y = y5;
				fifth_prediction[user_count][t].y = y6;

			}
			if (t == 2)
			{
				first_prediction[user_count][t].y = y3;
				second_prediction[user_count][t].y = y4;
				third_prediction[user_count][t].y = y5;
				fourth_prediction[user_count][t].y = y6;
				fifth_prediction[user_count][t].y = y7;

			}
			if (t == 3)
			{
				first_prediction[user_count][t].y = y4;
				second_prediction[user_count][t].y = y5;
				third_prediction[user_count][t].y = y6;
				fourth_prediction[user_count][t].y = y7;
				fifth_prediction[user_count][t].y = y8;

			}
			if (t == 4)
			{
				first_prediction[user_count][t].y = y5;
				second_prediction[user_count][t].y = y6;
				third_prediction[user_count][t].y = y7;
				fourth_prediction[user_count][t].y = y8;
				fifth_prediction[user_count][t].y = y9;

			}
			if (t == 5)
			{
				first_prediction[user_count][t].y = y6;
				second_prediction[user_count][t].y = y7;
				third_prediction[user_count][t].y = y8;
				fourth_prediction[user_count][t].y = y9;
				fifth_prediction[user_count][t].y = y10;

			}
			if (t == 6)
			{
				first_prediction[user_count][t].y = y7;
				second_prediction[user_count][t].y = y8;
				third_prediction[user_count][t].y = y9;
				fourth_prediction[user_count][t].y = y10;
				fifth_prediction[user_count][t].y = y11;

			}
			if (t == 7)
			{
				first_prediction[user_count][t].y = y8;
				second_prediction[user_count][t].y = y9;
				third_prediction[user_count][t].y = y10;
				fourth_prediction[user_count][t].y = y11;
				fifth_prediction[user_count][t].y = y12;

			}
			if (t == 8)
			{
				first_prediction[user_count][t].y = y9;
				second_prediction[user_count][t].y = y10;
				third_prediction[user_count][t].y = y11;
				fourth_prediction[user_count][t].y = y12;
				fifth_prediction[user_count][t].y = y13;

			}
			if (t == 9)
			{
				first_prediction[user_count][t].y = y10;
				second_prediction[user_count][t].y = y11;
				third_prediction[user_count][t].y = y12;
				fourth_prediction[user_count][t].y = y13;
				fifth_prediction[user_count][t].y = y14;

			}
			if (t == 10)
			{
				first_prediction[user_count][t].y = y11;
				second_prediction[user_count][t].y = y12;
				third_prediction[user_count][t].y = y13;
				fourth_prediction[user_count][t].y = y14;
				fifth_prediction[user_count][t].y = y15;

			}
			if (t == 11)
			{
				first_prediction[user_count][t].y = y12;
				second_prediction[user_count][t].y = y13;
				third_prediction[user_count][t].y = y14;
				fourth_prediction[user_count][t].y = y15;
				fifth_prediction[user_count][t].y = y16;

			}
			if (t == 12)
			{
				first_prediction[user_count][t].y = y13;
				second_prediction[user_count][t].y = y14;
				third_prediction[user_count][t].y = y15;
				fourth_prediction[user_count][t].y = y16;
				fifth_prediction[user_count][t].y = y17;

			}
			if (t == 13)
			{
				first_prediction[user_count][t].y = y14;
				second_prediction[user_count][t].y = y15;
				third_prediction[user_count][t].y = y16;
				fourth_prediction[user_count][t].y = y17;
				fifth_prediction[user_count][t].y = y18;

			}
			if (t == 14)
			{
				first_prediction[user_count][t].y = y15;
				second_prediction[user_count][t].y = y16;
				third_prediction[user_count][t].y = y17;
				fourth_prediction[user_count][t].y = y18;
				fifth_prediction[user_count][t].y = y19;

			}
			if (t == 15)
			{
				first_prediction[user_count][t].y = y16;
				second_prediction[user_count][t].y = y17;
				third_prediction[user_count][t].y = y18;
				fourth_prediction[user_count][t].y = y19;
				fifth_prediction[user_count][t].y = y20;

			}
			if (t == 16)
			{
				first_prediction[user_count][t].y = y17;
				second_prediction[user_count][t].y = y18;
				third_prediction[user_count][t].y = y19;
				fourth_prediction[user_count][t].y = y20;
			}
			if (t == 17)
			{
				first_prediction[user_count][t].y = y18;
				second_prediction[user_count][t].y = y19;
				third_prediction[user_count][t].y = y20;

			}
			if (t == 18)
			{
				first_prediction[user_count][t].y = y19;
				second_prediction[user_count][t].y = y20;

			}
			if (t == 19)
			{
				first_prediction[user_count][t].y = y20;

			}


			user_count++;
		}
		t++;
	}
	}
	/* id */
	{
	int t = 1;
	while (t <= 19)
	{
		stringstream out;
		out << t;
		
		ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\id\\Results_50\\results_Sample_Data_50_id_" + out.str() + ".csv");

		if (!data.is_open()) cout << "ERROR: File Open" << '\n';

		string one;
		string two;
		string three;
		string four;
		string five;
		string six;
		string seven;
		string eight;
		string nine;
		string ten;
		string eleven;
		string twelve;
		string thirteen;
		string fourteen;
		string fifteen;
		string sixteen;
		string seventeen;
		string eighteen;
		string nineteen;
		string twenty;

		int user_count = 0;
		
		while (data.good() && user_count<user_total)
		{

			getline(data, one, ',');
			getline(data, two, ',');
			getline(data, three, ',');
			getline(data, four, ',');
			getline(data, five, ',');
			getline(data, six, ',');
			getline(data, seven, ',');
			getline(data, eight, ',');
			getline(data, nine, ',');
			getline(data, ten, ',');
			getline(data, eleven, ',');
			getline(data, twelve, ',');
			getline(data, thirteen, ',');
			getline(data, fourteen, ',');
			getline(data, fifteen, ',');
			getline(data, sixteen, ',');
			getline(data, seventeen, ',');
			getline(data, eighteen, ',');
			getline(data, nineteen, ',');
			getline(data, twenty, '\n');



			/*  ---------------------------   */
			stringstream geek(one);
			double id1 = 0;
			geek >> id1;

			stringstream geek2(two);
			double id2 = 0;
			geek2 >> id2;

			stringstream geek3(three);
			double id3 = 0;
			geek3 >> id3;

			stringstream geek4(four);
			double id4 = 0;
			geek4 >> id4;

			stringstream geek5(five);
			double id5 = 0;
			geek5 >> id5;


			stringstream geek6(six);
			double id6 = 0;
			geek6 >> id6;

			stringstream geek7(seven);
			double id7 = 0;
			geek7 >> id7;

			stringstream geek8(eight);
			double id8 = 0;
			geek8 >> id8;

			stringstream geek9(nine);
			double id9 = 0;
			geek9 >> id9;

			stringstream geek10(ten);
			double id10 = 0;
			geek10 >> id10;

			stringstream geek11(eleven);
			double id11 = 0;
			geek11 >> id11;

			stringstream geek12(twelve);
			double id12 = 0;
			geek12 >> id12;

			stringstream geek13(thirteen);
			double id13 = 0;
			geek13 >> id13;

			stringstream geek14(fourteen);
			double id14 = 0;
			geek14 >> id14;

			stringstream geek15(fifteen);
			double id15 = 0;
			geek15 >> id15;

			stringstream geek16(sixteen);
			double id16 = 0;
			geek16 >> id16;

			stringstream geek17(seventeen);
			double id17 = 0;
			geek17 >> id17;

			stringstream geek18(eighteen);
			double id18 = 0;
			geek18 >> id18;

			stringstream geek19(nineteen);
			double id19 = 0;
			geek19 >> id19;

			stringstream geek20(twenty);
			double id20 = 0;
			geek20 >> id20;

			if (t == 1)
			{
				first_prediction[user_count][t].id = id2;
				second_prediction[user_count][t].id = id3;
				third_prediction[user_count][t].id = id4;
				fourth_prediction[user_count][t].id = id5;
				fifth_prediction[user_count][t].id = id6;

			}
			if (t == 2)
			{
				first_prediction[user_count][t].id = id3;
				second_prediction[user_count][t].id = id4;
				third_prediction[user_count][t].id = id5;
				fourth_prediction[user_count][t].id = id6;
				fifth_prediction[user_count][t].id = id7;

			}
			if (t == 3)
			{
				first_prediction[user_count][t].id = id4;
				second_prediction[user_count][t].id = id5;
				third_prediction[user_count][t].id = id6;
				fourth_prediction[user_count][t].id = id7;
				fifth_prediction[user_count][t].id = id8;

			}
			if (t == 4)
			{
				first_prediction[user_count][t].id = id5;
				second_prediction[user_count][t].id = id6;
				third_prediction[user_count][t].id = id7;
				fourth_prediction[user_count][t].id = id8;
				fifth_prediction[user_count][t].id = id9;

			}
			if (t == 5)
			{
				first_prediction[user_count][t].id = id6;
				second_prediction[user_count][t].id = id7;
				third_prediction[user_count][t].id = id8;
				fourth_prediction[user_count][t].id = id9;
				fifth_prediction[user_count][t].id = id10;

			}
			if (t == 6)
			{
				first_prediction[user_count][t].id = id7;
				second_prediction[user_count][t].id = id8;
				third_prediction[user_count][t].id = id9;
				fourth_prediction[user_count][t].id = id10;
				fifth_prediction[user_count][t].id = id11;

			}
			if (t == 7)
			{
				first_prediction[user_count][t].id = id8;
				second_prediction[user_count][t].id = id9;
				third_prediction[user_count][t].id = id10;
				fourth_prediction[user_count][t].id = id11;
				fifth_prediction[user_count][t].id = id12;

			}
			if (t == 8)
			{
				first_prediction[user_count][t].id = id9;
				second_prediction[user_count][t].id = id10;
				third_prediction[user_count][t].id = id11;
				fourth_prediction[user_count][t].id = id12;
				fifth_prediction[user_count][t].id = id13;

			}
			if (t == 9)
			{
				first_prediction[user_count][t].id = id10;
				second_prediction[user_count][t].id = id11;
				third_prediction[user_count][t].id = id12;
				fourth_prediction[user_count][t].id = id13;
				fifth_prediction[user_count][t].id = id14;

			}
			if (t == 10)
			{
				first_prediction[user_count][t].id = id11;
				second_prediction[user_count][t].id = id12;
				third_prediction[user_count][t].id = id13;
				fourth_prediction[user_count][t].id = id14;
				fifth_prediction[user_count][t].id = id15;

			}
			if (t == 11)
			{
				first_prediction[user_count][t].id = id12;
				second_prediction[user_count][t].id = id13;
				third_prediction[user_count][t].id = id14;
				fourth_prediction[user_count][t].id = id15;
				fifth_prediction[user_count][t].id = id16;

			}
			if (t == 12)
			{
				first_prediction[user_count][t].id = id13;
				second_prediction[user_count][t].id = id14;
				third_prediction[user_count][t].id = id15;
				fourth_prediction[user_count][t].id = id16;
				fifth_prediction[user_count][t].id = id17;

			}
			if (t == 13)
			{
				first_prediction[user_count][t].id = id14;
				second_prediction[user_count][t].id = id15;
				third_prediction[user_count][t].id = id16;
				fourth_prediction[user_count][t].id = id17;
				fifth_prediction[user_count][t].id = id18;

			}
			if (t == 14)
			{
				first_prediction[user_count][t].id = id15;
				second_prediction[user_count][t].id = id16;
				third_prediction[user_count][t].id = id17;
				fourth_prediction[user_count][t].id = id18;
				fifth_prediction[user_count][t].id = id19;

			}
			if (t == 15)
			{
				first_prediction[user_count][t].id = id16;
				second_prediction[user_count][t].id = id17;
				third_prediction[user_count][t].id = id18;
				fourth_prediction[user_count][t].id = id19;
				fifth_prediction[user_count][t].id = id20;

			}
			if (t == 16)
			{
				first_prediction[user_count][t].id = id17;
				second_prediction[user_count][t].id = id18;
				third_prediction[user_count][t].id = id19;
				fourth_prediction[user_count][t].id = id20;

			}
			if (t == 17)
			{
				first_prediction[user_count][t].id = id18;
				second_prediction[user_count][t].id = id19;
				third_prediction[user_count][t].id = id20;

			}
			if (t == 18)
			{
				first_prediction[user_count][t].id = id19;
				second_prediction[user_count][t].id = id20;

			}
			if (t == 19)
			{
				first_prediction[user_count][t].id = id20;

			}


			user_count++;
		}
		t++;
	}
	}
	/* md */
	{
	int t = 1;
	while (t <= 19)
	{
		stringstream out;
		out << t;

		ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\md\\Results_50\\results_Sample_Data_50_md_" + out.str() + ".csv");

		if (!data.is_open()) cout << "ERROR: File Open" << '\n';

		string one;
		string two;
		string three;
		string four;
		string five;
		string six;
		string seven;
		string eight;
		string nine;
		string ten;
		string eleven;
		string twelve;
		string thirteen;
		string fourteen;
		string fifteen;
		string sixteen;
		string seventeen;
		string eighteen;
		string nineteen;
		string twenty;

		int user_count = 0;
	
		while (data.good() && user_count<user_total)
		{

			getline(data, one, ',');
			getline(data, two, ',');
			getline(data, three, ',');
			getline(data, four, ',');
			getline(data, five, ',');
			getline(data, six, ',');
			getline(data, seven, ',');
			getline(data, eight, ',');
			getline(data, nine, ',');
			getline(data, ten, ',');
			getline(data, eleven, ',');
			getline(data, twelve, ',');
			getline(data, thirteen, ',');
			getline(data, fourteen, ',');
			getline(data, fifteen, ',');
			getline(data, sixteen, ',');
			getline(data, seventeen, ',');
			getline(data, eighteen, ',');
			getline(data, nineteen, ',');
			getline(data, twenty, '\n');



			/*  ---------------------------   */
			stringstream geek(one);
			double md1 = 0;
			geek >> md1;

			stringstream geek2(two);
			double md2 = 0;
			geek2 >> md2;

			stringstream geek3(three);
			double md3 = 0;
			geek3 >> md3;

			stringstream geek4(four);
			double md4 = 0;
			geek4 >> md4;

			stringstream geek5(five);
			double md5 = 0;
			geek5 >> md5;

			stringstream geek6(six);
			double md6 = 0;
			geek6 >> md6;

			stringstream geek7(seven);
			double md7 = 0;
			geek7 >> md7;

			stringstream geek8(eight);
			double md8 = 0;
			geek8 >> md8;

			stringstream geek9(nine);
			double md9 = 0;
			geek9 >> md9;

			stringstream geek10(ten);
			double md10 = 0;
			geek10 >> md10;

			stringstream geek11(eleven);
			double md11 = 0;
			geek11 >> md11;

			stringstream geek12(twelve);
			double md12 = 0;
			geek12 >> md12;

			stringstream geek13(thirteen);
			double md13 = 0;
			geek13 >> md13;

			stringstream geek14(fourteen);
			double md14 = 0;
			geek14 >> md14;

			stringstream geek15(fifteen);
			double md15 = 0;
			geek15 >> md15;

			stringstream geek16(sixteen);
			double md16 = 0;
			geek16 >> md16;

			stringstream geek17(seventeen);
			double md17 = 0;
			geek17 >> md17;

			stringstream geek18(eighteen);
			double md18 = 0;
			geek18 >> md18;

			stringstream geek19(nineteen);
			double md19 = 0;
			geek19 >> md19;

			stringstream geek20(twenty);
			double md20 = 0;
			geek20 >> md20;

			if (t == 1)
			{
				first_prediction[user_count][t].md = md2;
				second_prediction[user_count][t].md = md3;
				third_prediction[user_count][t].md = md4;
				fourth_prediction[user_count][t].md = md5;
				fifth_prediction[user_count][t].md = md6;

			}
			if (t == 2)
			{
				first_prediction[user_count][t].md = md3;
				second_prediction[user_count][t].md = md4;
				third_prediction[user_count][t].md = md5;
				fourth_prediction[user_count][t].md = md6;
				fifth_prediction[user_count][t].md = md7;

			}
			if (t == 3)
			{
				first_prediction[user_count][t].md = md4;
				second_prediction[user_count][t].md = md5;
				third_prediction[user_count][t].md = md6;
				fourth_prediction[user_count][t].md = md7;
				fifth_prediction[user_count][t].md = md8;

			}
			if (t == 4)
			{
				first_prediction[user_count][t].md = md5;
				second_prediction[user_count][t].md = md6;
				third_prediction[user_count][t].md = md7;
				fourth_prediction[user_count][t].md = md8;
				fifth_prediction[user_count][t].md = md9;

			}
			if (t == 5)
			{
				first_prediction[user_count][t].md = md6;
				second_prediction[user_count][t].md = md7;
				third_prediction[user_count][t].md = md8;
				fourth_prediction[user_count][t].md = md9;
				fifth_prediction[user_count][t].md = md10;

			}
			if (t == 6)
			{
				first_prediction[user_count][t].md = md7;
				second_prediction[user_count][t].md = md8;
				third_prediction[user_count][t].md = md9;
				fourth_prediction[user_count][t].md = md10;
				fifth_prediction[user_count][t].md = md11;

			}
			if (t == 7)
			{
				first_prediction[user_count][t].md = md8;
				second_prediction[user_count][t].md = md9;
				third_prediction[user_count][t].md = md10;
				fourth_prediction[user_count][t].md = md11;
				fifth_prediction[user_count][t].md = md12;

			}
			if (t == 8)
			{
				first_prediction[user_count][t].md = md9;
				second_prediction[user_count][t].md = md10;
				third_prediction[user_count][t].md = md11;
				fourth_prediction[user_count][t].md = md12;
				fifth_prediction[user_count][t].md = md13;

			}
			if (t == 9)
			{
				first_prediction[user_count][t].md = md10;
				second_prediction[user_count][t].md = md11;
				third_prediction[user_count][t].md = md12;
				fourth_prediction[user_count][t].md = md13;
				fifth_prediction[user_count][t].md = md14;

			}
			if (t == 10)
			{
				first_prediction[user_count][t].md = md11;
				second_prediction[user_count][t].md = md12;
				third_prediction[user_count][t].md = md13;
				fourth_prediction[user_count][t].md = md14;
				fifth_prediction[user_count][t].md = md15;

			}
			if (t == 11)
			{
				first_prediction[user_count][t].md = md12;
				second_prediction[user_count][t].md = md13;
				third_prediction[user_count][t].md = md14;
				fourth_prediction[user_count][t].md = md15;
				fifth_prediction[user_count][t].md = md16;

			}
			if (t == 12)
			{
				first_prediction[user_count][t].md = md13;
				second_prediction[user_count][t].md = md14;
				third_prediction[user_count][t].md = md15;
				fourth_prediction[user_count][t].md = md16;
				fifth_prediction[user_count][t].md = md17;

			}
			if (t == 13)
			{
				first_prediction[user_count][t].md = md14;
				second_prediction[user_count][t].md = md15;
				third_prediction[user_count][t].md = md16;
				fourth_prediction[user_count][t].md = md17;
				fifth_prediction[user_count][t].md = md18;

			}
			if (t == 14)
			{
				first_prediction[user_count][t].md = md15;
				second_prediction[user_count][t].md = md16;
				third_prediction[user_count][t].md = md17;
				fourth_prediction[user_count][t].md = md18;
				fifth_prediction[user_count][t].md = md19;

			}
			if (t == 15)
			{
				first_prediction[user_count][t].md = md16;
				second_prediction[user_count][t].md = md17;
				third_prediction[user_count][t].md = md18;
				fourth_prediction[user_count][t].md = md19;
				fifth_prediction[user_count][t].md = md20;

			}
			if (t == 16)
			{
				first_prediction[user_count][t].md = md17;
				second_prediction[user_count][t].md = md18;
				third_prediction[user_count][t].md = md19;
				fourth_prediction[user_count][t].md = md20;

			}
			if (t == 17)
			{
				first_prediction[user_count][t].md = md18;
				second_prediction[user_count][t].md = md19;
				third_prediction[user_count][t].md = md20;
			}
			if (t == 18)
			{
				first_prediction[user_count][t].md = md19;
				second_prediction[user_count][t].md = md20;

			}
			if (t == 19)
			{
				first_prediction[user_count][t].md = md20;

			}


			user_count++;
		}
		t++;
	}
	}
   /* w */
   {
	int t = 1;
	while (t <= 19)
	{
		stringstream out;
		out << t;

		ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\w\\Results_50\\results_Sample_Data_50_w_" + out.str() + ".csv");

		if (!data.is_open()) cout << "ERROR: File Open" << '\n';

		string one;
		string two;
		string three;
		string four;
		string five;
		string six;
		string seven;
		string eight;
		string nine;
		string ten;
		string eleven;
		string twelve;
		string thirteen;
		string fourteen;
		string fifteen;
		string sixteen;
		string seventeen;
		string eighteen;
		string nineteen;
		string twenty;

		int user_count = 0;
		
		while (data.good() && user_count<user_total)
		{

			getline(data, one, ',');
			getline(data, two, ',');
			getline(data, three, ',');
			getline(data, four, ',');
			getline(data, five, ',');
			getline(data, six, ',');
			getline(data, seven, ',');
			getline(data, eight, ',');
			getline(data, nine, ',');
			getline(data, ten, ',');
			getline(data, eleven, ',');
			getline(data, twelve, ',');
			getline(data, thirteen, ',');
			getline(data, fourteen, ',');
			getline(data, fifteen, ',');
			getline(data, sixteen, ',');
			getline(data, seventeen, ',');
			getline(data, eighteen, ',');
			getline(data, nineteen, ',');
			getline(data, twenty, '\n');



			/*  ---------------------------   */
			stringstream geek(one);
			double w1 = 0;
			geek >> w1;

			stringstream geek2(two);
			double w2 = 0;
			geek2 >> w2;

			stringstream geek3(three);
			double w3 = 0;
			geek3 >> w3;

			stringstream geek4(four);
			double w4 = 0;
			geek4 >> w4;

			stringstream geek5(five);
			double w5 = 0;
			geek5 >> w5;

			stringstream geek6(six);
			double w6 = 0;
			geek6 >> w6;

			stringstream geek7(seven);
			double w7 = 0;
			geek7 >> w7;

			stringstream geek8(eight);
			double w8 = 0;
			geek8 >> w8;

			stringstream geek9(nine);
			double w9 = 0;
			geek9 >> w9;

			stringstream geek10(ten);
			double w10 = 0;
			geek10 >> w10;

			stringstream geek11(eleven);
			double w11 = 0;
			geek11 >> w11;

			stringstream geek12(twelve);
			double w12 = 0;
			geek12 >> w12;

			stringstream geek13(thirteen);
			double w13 = 0;
			geek13 >> w13;

			stringstream geek14(fourteen);
			double w14 = 0;
			geek14 >> w14;

			stringstream geek15(fifteen);
			double w15 = 0;
			geek15 >> w15;

			stringstream geek16(sixteen);
			double w16 = 0;
			geek16 >> w16;

			stringstream geek17(seventeen);
			double w17 = 0;
			geek17 >> w17;

			stringstream geek18(eighteen);
			double w18 = 0;
			geek18 >> w18;

			stringstream geek19(nineteen);
			double w19 = 0;
			geek19 >> w19;

			stringstream geek20(twenty);
			double w20 = 0;
			geek20 >> w20;

			if (t == 1)
			{
				first_prediction[user_count][t].w = w2;
				second_prediction[user_count][t].w = w3;
				third_prediction[user_count][t].w = w4;
				fourth_prediction[user_count][t].w = w5;
				fifth_prediction[user_count][t].w = w6;

			}
			if (t == 2)
			{
				first_prediction[user_count][t].w = w3;
				second_prediction[user_count][t].w = w4;
				third_prediction[user_count][t].w = w5;
				fourth_prediction[user_count][t].w = w6;
				fifth_prediction[user_count][t].w = w7;

			}
			if (t == 3)
			{
				first_prediction[user_count][t].w = w4;
				second_prediction[user_count][t].w = w5;
				third_prediction[user_count][t].w = w6;
				fourth_prediction[user_count][t].w = w7;
				fifth_prediction[user_count][t].w = w8;

			}
			if (t == 4)
			{
				first_prediction[user_count][t].w = w5;
				second_prediction[user_count][t].w = w6;
				third_prediction[user_count][t].w = w7;
				fourth_prediction[user_count][t].w = w8;
				fifth_prediction[user_count][t].w = w9;

			}
			if (t == 5)
			{
				first_prediction[user_count][t].w = w6;
				second_prediction[user_count][t].w = w7;
				third_prediction[user_count][t].w = w8;
				fourth_prediction[user_count][t].w = w9;
				fifth_prediction[user_count][t].w = w10;

			}
			if (t == 6)
			{
				first_prediction[user_count][t].w = w7;
				second_prediction[user_count][t].w = w8;
				third_prediction[user_count][t].w = w9;
				fourth_prediction[user_count][t].w = w10;
				fifth_prediction[user_count][t].w = w11;

			}
			if (t == 7)
			{
				first_prediction[user_count][t].w = w8;
				second_prediction[user_count][t].w = w9;
				third_prediction[user_count][t].w = w10;
				fourth_prediction[user_count][t].w = w11;
				fifth_prediction[user_count][t].w = w12;

			}
			if (t == 8)
			{
				first_prediction[user_count][t].w = w9;
				second_prediction[user_count][t].w = w10;
				third_prediction[user_count][t].w = w11;
				fourth_prediction[user_count][t].w = w12;
				fifth_prediction[user_count][t].w = w13;

			}
			if (t == 9)
			{
				first_prediction[user_count][t].w = w10;
				second_prediction[user_count][t].w = w11;
				third_prediction[user_count][t].w = w12;
				fourth_prediction[user_count][t].w = w13;
				fifth_prediction[user_count][t].w = w14;

			}
			if (t == 10)
			{
				first_prediction[user_count][t].w = w11;
				second_prediction[user_count][t].w = w12;
				third_prediction[user_count][t].w = w13;
				fourth_prediction[user_count][t].w = w14;
				fifth_prediction[user_count][t].w = w15;

			}
			if (t == 11)
			{
				first_prediction[user_count][t].w = w12;
				second_prediction[user_count][t].w = w13;
				third_prediction[user_count][t].w = w14;
				fourth_prediction[user_count][t].w = w15;
				fifth_prediction[user_count][t].w = w16;

			}
			if (t == 12)
			{
				first_prediction[user_count][t].w = w13;
				second_prediction[user_count][t].w = w14;
				third_prediction[user_count][t].w = w15;
				fourth_prediction[user_count][t].w = w16;
				fifth_prediction[user_count][t].w = w17;

			}
			if (t == 13)
			{
				first_prediction[user_count][t].w = w14;
				second_prediction[user_count][t].w = w15;
				third_prediction[user_count][t].w = w16;
				fourth_prediction[user_count][t].w = w17;
				fifth_prediction[user_count][t].w = w18;

			}
			if (t == 14)
			{
				first_prediction[user_count][t].w = w15;
				second_prediction[user_count][t].w = w16;
				third_prediction[user_count][t].w = w17;
				fourth_prediction[user_count][t].w = w18;
				fifth_prediction[user_count][t].w = w19;

			}
			if (t == 15)
			{
				first_prediction[user_count][t].w = w16;
				second_prediction[user_count][t].w = w17;
				third_prediction[user_count][t].w = w18;
				fourth_prediction[user_count][t].w = w19;
				fifth_prediction[user_count][t].w = w20;

			}
			if (t == 16)
			{
				first_prediction[user_count][t].w = w17;
				second_prediction[user_count][t].w = w18;
				third_prediction[user_count][t].w = w19;
				fourth_prediction[user_count][t].w = w20;
			}
			if (t == 17)
			{
				first_prediction[user_count][t].w = w18;
				second_prediction[user_count][t].w = w19;
				third_prediction[user_count][t].w = w20;
			}
			if (t == 18)
			{
				first_prediction[user_count][t].w = w19;
				second_prediction[user_count][t].w = w20;
			}
			if (t == 19)
			{
				first_prediction[user_count][t].w = w20;

			}


			user_count++;
		}
		t++;
	}
	}
   /* p */
   {
	int t = 1;
	while (t <= 19)
	{
		stringstream out;
		out << t;
		
		ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\p\\Results_50\\results_Sample_Data_50_p_" + out.str() + ".csv");

		if (!data.is_open()) cout << "ERROR: File Open" << '\n';

		string one;
		string two;
		string three;
		string four;
		string five;
		string six;
		string seven;
		string eight;
		string nine;
		string ten;
		string eleven;
		string twelve;
		string thirteen;
		string fourteen;
		string fifteen;
		string sixteen;
		string seventeen;
		string eighteen;
		string nineteen;
		string twenty;

		int user_count = 0;
	
		while (data.good() && user_count<user_total)
		{

			getline(data, one, ',');
			getline(data, two, ',');
			getline(data, three, ',');
			getline(data, four, ',');
			getline(data, five, ',');
			getline(data, six, ',');
			getline(data, seven, ',');
			getline(data, eight, ',');
			getline(data, nine, ',');
			getline(data, ten, ',');
			getline(data, eleven, ',');
			getline(data, twelve, ',');
			getline(data, thirteen, ',');
			getline(data, fourteen, ',');
			getline(data, fifteen, ',');
			getline(data, sixteen, ',');
			getline(data, seventeen, ',');
			getline(data, eighteen, ',');
			getline(data, nineteen, ',');
			getline(data, twenty, '\n');



			/*  ---------------------------   */
			stringstream geek(one);
			double p1 = 0;
			geek >> p1;

			stringstream geek2(two);
			double p2 = 0;
			geek2 >> p2;

			stringstream geek3(three);
			double p3 = 0;
			geek3 >> p3;

			stringstream geek4(four);
			double p4 = 0;
			geek4 >> p4;

			stringstream geek5(five);
			double p5 = 0;
			geek5 >> p5;

			stringstream geek6(six);
			double p6 = 0;
			geek6 >> p6;

			stringstream geek7(seven);
			double p7 = 0;
			geek7 >> p7;

			stringstream geek8(eight);
			double p8 = 0;
			geek8 >> p8;

			stringstream geek9(nine);
			double p9 = 0;
			geek9 >> p9;

			stringstream geek10(ten);
			double p10 = 0;
			geek10 >> p10;

			stringstream geek11(eleven);
			double p11 = 0;
			geek11 >> p11;

			stringstream geek12(twelve);
			double p12 = 0;
			geek12 >> p12;

			stringstream geek13(thirteen);
			double p13 = 0;
			geek13 >> p13;

			stringstream geek14(fourteen);
			double p14 = 0;
			geek14 >> p14;

			stringstream geek15(fifteen);
			double p15 = 0;
			geek15 >> p15;

			stringstream geek16(sixteen);
			double p16 = 0;
			geek16 >> p16;

			stringstream geek17(seventeen);
			double p17 = 0;
			geek17 >> p17;

			stringstream geek18(eighteen);
			double p18 = 0;
			geek18 >> p18;

			stringstream geek19(nineteen);
			double p19 = 0;
			geek19 >> p19;

			stringstream geek20(twenty);
			double p20 = 0;
			geek20 >> p20;

			if (t == 1)
			{
				first_prediction[user_count][t].p = p2;
				second_prediction[user_count][t].p = p3;
				third_prediction[user_count][t].p = p4;
				fourth_prediction[user_count][t].p = p5;
				fifth_prediction[user_count][t].p = p6;

			}
			if (t == 2)
			{
				first_prediction[user_count][t].p = p3;
				second_prediction[user_count][t].p = p4;
				third_prediction[user_count][t].p = p5;
				fourth_prediction[user_count][t].p = p6;
				fifth_prediction[user_count][t].p = p7;

			}
			if (t == 3)
			{
				first_prediction[user_count][t].p = p4;
				second_prediction[user_count][t].p = p5;
				third_prediction[user_count][t].p = p6;
				fourth_prediction[user_count][t].p = p7;
				fifth_prediction[user_count][t].p = p8;

			}
			if (t == 4)
			{
				first_prediction[user_count][t].p = p5;
				second_prediction[user_count][t].p = p6;
				third_prediction[user_count][t].p = p7;
				fourth_prediction[user_count][t].p = p8;
				fifth_prediction[user_count][t].p = p9;

			}
			if (t == 5)
			{
				first_prediction[user_count][t].p = p6;
				second_prediction[user_count][t].p = p7;
				third_prediction[user_count][t].p = p8;
				fourth_prediction[user_count][t].p = p9;
				fifth_prediction[user_count][t].p = p10;

			}
			if (t == 6)
			{
				first_prediction[user_count][t].p = p7;
				second_prediction[user_count][t].p = p8;
				third_prediction[user_count][t].p = p9;
				fourth_prediction[user_count][t].p = p10;
				fifth_prediction[user_count][t].p = p11;

			}
			if (t == 7)
			{
				first_prediction[user_count][t].p = p8;
				second_prediction[user_count][t].p = p9;
				third_prediction[user_count][t].p = p10;
				fourth_prediction[user_count][t].p = p11;
				fifth_prediction[user_count][t].p = p12;

			}
			if (t == 8)
			{
				first_prediction[user_count][t].p = p9;
				second_prediction[user_count][t].p = p10;
				third_prediction[user_count][t].p = p11;
				fourth_prediction[user_count][t].p = p12;
				fifth_prediction[user_count][t].p = p13;

			}
			if (t == 9)
			{
				first_prediction[user_count][t].p = p10;
				second_prediction[user_count][t].p = p11;
				third_prediction[user_count][t].p = p12;
				fourth_prediction[user_count][t].p = p13;
				fifth_prediction[user_count][t].p = p14;

			}
			if (t == 10)
			{
				first_prediction[user_count][t].p = p11;
				second_prediction[user_count][t].p = p12;
				third_prediction[user_count][t].p = p13;
				fourth_prediction[user_count][t].p = p14;
				fifth_prediction[user_count][t].p = p15;

			}
			if (t == 11)
			{
				first_prediction[user_count][t].p = p12;
				second_prediction[user_count][t].p = p13;
				third_prediction[user_count][t].p = p14;
				fourth_prediction[user_count][t].p = p15;
				fifth_prediction[user_count][t].p = p16;

			}
			if (t == 12)
			{
				first_prediction[user_count][t].p = p13;
				second_prediction[user_count][t].p = p14;
				third_prediction[user_count][t].p = p15;
				fourth_prediction[user_count][t].p = p16;
				fifth_prediction[user_count][t].p = p17;

			}
			if (t == 13)
			{
				first_prediction[user_count][t].p = p14;
				second_prediction[user_count][t].p = p15;
				third_prediction[user_count][t].p = p16;
				fourth_prediction[user_count][t].p = p17;
				fifth_prediction[user_count][t].p = p18;

			}
			if (t == 14)
			{
				first_prediction[user_count][t].p = p15;
				second_prediction[user_count][t].p = p16;
				third_prediction[user_count][t].p = p17;
				fourth_prediction[user_count][t].p = p18;
				fifth_prediction[user_count][t].p = p19;

			}
			if (t == 15)
			{
				first_prediction[user_count][t].p = p16;
				second_prediction[user_count][t].p = p17;
				third_prediction[user_count][t].p = p18;
				fourth_prediction[user_count][t].p = p19;
				fifth_prediction[user_count][t].p = p20;

			}
			if (t == 16)
			{
				first_prediction[user_count][t].p = p17;
				second_prediction[user_count][t].p = p18;
				third_prediction[user_count][t].p = p19;
				fourth_prediction[user_count][t].p = p20;

			}
			if (t == 17)
			{
				first_prediction[user_count][t].p = p18;
				second_prediction[user_count][t].p = p19;
				third_prediction[user_count][t].p = p20;
			}
			if (t == 18)
			{
				first_prediction[user_count][t].p = p19;
				second_prediction[user_count][t].p = p20;
			}
			if (t == 19)
			{
				first_prediction[user_count][t].p = p20;

			}


			user_count++;
		}
		t++;
	}
    }
	/* b */
	{
	int t = 1;
	while (t <= 19)
	{
		stringstream out;
		out << t;
		
		ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\b\\Results_50\\results_Sample_Data_50_b_" + out.str() + ".csv");

		if (!data.is_open()) cout << "ERROR: File Open" << '\n';

		string one;
		string two;
		string three;
		string four;
		string five;
		string six;
		string seven;
		string eight;
		string nine;
		string ten;
		string eleven;
		string twelve;
		string thirteen;
		string fourteen;
		string fifteen;
		string sixteen;
		string seventeen;
		string eighteen;
		string nineteen;
		string twenty;

		int user_count = 0;
		
		while (data.good() && user_count<user_total)
		{

			getline(data, one, ',');
			getline(data, two, ',');
			getline(data, three, ',');
			getline(data, four, ',');
			getline(data, five, ',');
			getline(data, six, ',');
			getline(data, seven, ',');
			getline(data, eight, ',');
			getline(data, nine, ',');
			getline(data, ten, ',');
			getline(data, eleven, ',');
			getline(data, twelve, ',');
			getline(data, thirteen, ',');
			getline(data, fourteen, ',');
			getline(data, fifteen, ',');
			getline(data, sixteen, ',');
			getline(data, seventeen, ',');
			getline(data, eighteen, ',');
			getline(data, nineteen, ',');
			getline(data, twenty, '\n');



			/*  ---------------------------   */
			stringstream geek(one);
			double b1 = 0;
			geek >> b1;

			stringstream geek2(two);
			double b2 = 0;
			geek2 >> b2;

			stringstream geek3(three);
			double b3 = 0;
			geek3 >> b3;

			stringstream geek4(four);
			double b4 = 0;
			geek4 >> b4;

			stringstream geek5(five);
			double b5 = 0;
			geek5 >> b5;

			stringstream geek6(six);
			double b6 = 0;
			geek6 >> b6;

			stringstream geek7(seven);
			double b7 = 0;
			geek7 >> b7;

			stringstream geek8(eight);
			double b8 = 0;
			geek8 >> b8;

			stringstream geek9(nine);
			double b9 = 0;
			geek9 >> b9;

			stringstream geek10(ten);
			double b10 = 0;
			geek10 >> b10;

			stringstream geek11(eleven);
			double b11 = 0;
			geek11 >> b11;

			stringstream geek12(twelve);
			double b12 = 0;
			geek12 >> b12;

			stringstream geek13(thirteen);
			double b13 = 0;
			geek13 >> b13;

			stringstream geek14(fourteen);
			double b14 = 0;
			geek14 >> b14;

			stringstream geek15(fifteen);
			double b15 = 0;
			geek15 >> b15;

			stringstream geek16(sixteen);
			double b16 = 0;
			geek16 >> b16;

			stringstream geek17(seventeen);
			double b17 = 0;
			geek17 >> b17;

			stringstream geek18(eighteen);
			double b18 = 0;
			geek18 >> b18;

			stringstream geek19(nineteen);
			double b19 = 0;
			geek19 >> b19;

			stringstream geek20(twenty);
			double b20 = 0;
			geek20 >> b20;

			if (t == 1)
			{
				first_prediction[user_count][t].b = b2;
				second_prediction[user_count][t].b = b3;
				third_prediction[user_count][t].b = b4;
				fourth_prediction[user_count][t].b = b5;
				fifth_prediction[user_count][t].b = b6;

			}
			if (t == 2)
			{
				first_prediction[user_count][t].b = b3;
				second_prediction[user_count][t].b = b4;
				third_prediction[user_count][t].b = b5;
				fourth_prediction[user_count][t].b = b6;
				fifth_prediction[user_count][t].b = b7;

			}
			if (t == 3)
			{
				first_prediction[user_count][t].b = b4;
				second_prediction[user_count][t].b = b5;
				third_prediction[user_count][t].b = b6;
				fourth_prediction[user_count][t].b = b7;
				fifth_prediction[user_count][t].b = b8;

			}
			if (t == 4)
			{
				first_prediction[user_count][t].b = b5;
				second_prediction[user_count][t].b = b6;
				third_prediction[user_count][t].b = b7;
				fourth_prediction[user_count][t].b = b8;
				fifth_prediction[user_count][t].b = b9;

			}
			if (t == 5)
			{
				first_prediction[user_count][t].b = b6;
				second_prediction[user_count][t].b = b7;
				third_prediction[user_count][t].b = b8;
				fourth_prediction[user_count][t].b = b9;
				fifth_prediction[user_count][t].b = b10;

			}
			if (t == 6)
			{
				first_prediction[user_count][t].b = b7;
				second_prediction[user_count][t].b = b8;
				third_prediction[user_count][t].b = b9;
				fourth_prediction[user_count][t].b = b10;
				fifth_prediction[user_count][t].b = b11;

			}
			if (t == 7)
			{
				first_prediction[user_count][t].b = b8;
				second_prediction[user_count][t].b = b9;
				third_prediction[user_count][t].b = b10;
				fourth_prediction[user_count][t].b = b11;
				fifth_prediction[user_count][t].b = b12;

			}
			if (t == 8)
			{
				first_prediction[user_count][t].b = b9;
				second_prediction[user_count][t].b = b10;
				third_prediction[user_count][t].b = b11;
				fourth_prediction[user_count][t].b = b12;
				fifth_prediction[user_count][t].b = b13;

			}
			if (t == 9)
			{
				first_prediction[user_count][t].b = b10;
				second_prediction[user_count][t].b = b11;
				third_prediction[user_count][t].b = b12;
				fourth_prediction[user_count][t].b = b13;
				fifth_prediction[user_count][t].b = b14;

			}
			if (t == 10)
			{
				first_prediction[user_count][t].b = b11;
				second_prediction[user_count][t].b = b12;
				third_prediction[user_count][t].b = b13;
				fourth_prediction[user_count][t].b = b14;
				fifth_prediction[user_count][t].b = b15;

			}
			if (t == 11)
			{
				first_prediction[user_count][t].b = b12;
				second_prediction[user_count][t].b = b13;
				third_prediction[user_count][t].b = b14;
				fourth_prediction[user_count][t].b = b15;
				fifth_prediction[user_count][t].b = b16;

			}
			if (t == 12)
			{
				first_prediction[user_count][t].b = b13;
				second_prediction[user_count][t].b = b14;
				third_prediction[user_count][t].b = b15;
				fourth_prediction[user_count][t].b = b16;
				fifth_prediction[user_count][t].b = b17;

			}
			if (t == 13)
			{
				first_prediction[user_count][t].b = b14;
				second_prediction[user_count][t].b = b15;
				third_prediction[user_count][t].b = b16;
				fourth_prediction[user_count][t].b = b17;
				fifth_prediction[user_count][t].b = b18;

			}
			if (t == 14)
			{
				first_prediction[user_count][t].b = b15;
				second_prediction[user_count][t].b = b16;
				third_prediction[user_count][t].b = b17;
				fourth_prediction[user_count][t].b = b18;
				fifth_prediction[user_count][t].b = b19;

			}
			if (t == 15)
			{
				first_prediction[user_count][t].b = b16;
				second_prediction[user_count][t].b = b17;
				third_prediction[user_count][t].b = b18;
				fourth_prediction[user_count][t].b = b19;
				fifth_prediction[user_count][t].b = b20;

			}
			if (t == 16)
			{
				first_prediction[user_count][t].b = b17;
				second_prediction[user_count][t].b = b18;
				third_prediction[user_count][t].b = b19;
				fourth_prediction[user_count][t].b = b20;

			}
			if (t == 17)
			{
				first_prediction[user_count][t].b = b18;
				second_prediction[user_count][t].b = b19;
				third_prediction[user_count][t].b = b20;

			}
			if (t == 18)
			{
				first_prediction[user_count][t].b = b19;
				second_prediction[user_count][t].b = b20;
			}
			if (t == 19)
			{
				first_prediction[user_count][t].b = b20;

			}


			user_count++;
		}
		t++;
	}
	}
	 
}

void ReadBig()
{

	ifstream data("C:\\Users\\Erfan\\Desktop\\vm_cpu_readings-file-1-of-195.csv");


	if (!data.is_open()) cout << "ERROR: File Open" << '\n';



	int t = 1;
	stringstream out1;
	out1 << t;
	ofstream myfile1;
	myfile1.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out1.str() + ".csv");
	myfile1.precision(25);

	t++;
	stringstream out2;
	out2 << t;
	ofstream myfile2;
	myfile2.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out2.str() + ".csv");
	myfile2.precision(25);

	t++;
	stringstream out3;
	out3 << t;
	ofstream myfile3;
	myfile3.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out3.str() + ".csv");
	myfile3.precision(25);

	t++;
	stringstream out4;
	out4 << t;
	ofstream myfile4;
	myfile4.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out4.str() + ".csv");
	myfile4.precision(25);

	t++;
	stringstream out5;
	out5 << t;
	ofstream myfile5;
	myfile5.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out5.str() + ".csv");
	myfile5.precision(25);

	t++;
	stringstream out6;
	out6 << t;
	ofstream myfile6;
	myfile6.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out6.str() + ".csv");
	myfile6.precision(25);

	t++;
	stringstream out7;
	out7 << t;
	ofstream myfile7;
	myfile7.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out7.str() + ".csv");
	myfile7.precision(25);

	t++;
	stringstream out8;
	out8 << t;
	ofstream myfile8;
	myfile8.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out8.str() + ".csv");
	myfile8.precision(25);

	t++;
	stringstream out9;
	out9 << t;
	ofstream myfile9;
	myfile9.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out9.str() + ".csv");
	myfile9.precision(25);

	t++;
	stringstream out10;
	out10 << t;
	ofstream myfile10;
	myfile10.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out10.str() + ".csv");
	myfile10.precision(25);

	t++;
	stringstream out11;
	out11 << t;
	ofstream myfile11;
	myfile11.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out11.str() + ".csv");
	myfile11.precision(25);

	t++;
	stringstream out12;
	out12 << t;
	ofstream myfile12;
	myfile12.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out12.str() + ".csv");
	myfile12.precision(25);

	t++;
	stringstream out13;
	out13 << t;
	ofstream myfile13;
	myfile13.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out13.str() + ".csv");
	myfile13.precision(25);

	t++;
	stringstream out14;
	out14 << t;
	ofstream myfile14;
	myfile14.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out14.str() + ".csv");
	myfile14.precision(25);

	t++;
	stringstream out15;
	out15 << t;
	ofstream myfile15;
	myfile15.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out15.str() + ".csv");
	myfile15.precision(25);

	t++;
	stringstream out16;
	out16 << t;
	ofstream myfile16;
	myfile16.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out16.str() + ".csv");
	myfile16.precision(25);

	t++;
	stringstream out17;
	out17 << t;
	ofstream myfile17;
	myfile17.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out17.str() + ".csv");
	myfile17.precision(25);

	t++;
	stringstream out18;
	out18 << t;
	ofstream myfile18;
	myfile18.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out18.str() + ".csv");
	myfile18.precision(25);

	t++;
	stringstream out19;
	out19 << t;
	ofstream myfile19;
	myfile19.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out19.str() + ".csv");
	myfile19.precision(25);

	t++;
	stringstream out20;
	out20 << t;
	ofstream myfile20;
	myfile20.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out20.str() + ".csv");
	myfile20.precision(25);

	t++;
	stringstream out21;
	out21 << t;
	ofstream myfile21;
	myfile21.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out21.str() + ".csv");
	myfile21.precision(25);

	t++;
	stringstream out22;
	out22 << t;
	ofstream myfile22;
	myfile22.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out22.str() + ".csv");
	myfile22.precision(25);

	t++;
	stringstream out23;
	out23 << t;
	ofstream myfile23;
	myfile23.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out23.str() + ".csv");
	myfile23.precision(25);

	t++;
	stringstream out24;
	out24 << t;
	ofstream myfile24;
	myfile24.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out24.str() + ".csv");
	myfile24.precision(25);

	t++;
	stringstream out25;
	out25 << t;
	ofstream myfile25;
	myfile25.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out25.str() + ".csv");
	myfile25.precision(25);

	t++;
	stringstream out26;
	out26 << t;
	ofstream myfile26;
	myfile26.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out26.str() + ".csv");
	myfile26.precision(25);

	t++;
	stringstream out27;
	out27 << t;
	ofstream myfile27;
	myfile27.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out27.str() + ".csv");
	myfile27.precision(25);

	t++;
	stringstream out28;
	out28 << t;
	ofstream myfile28;
	myfile28.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out28.str() + ".csv");
	myfile28.precision(25);

	t++;
	stringstream out29;
	out29 << t;
	ofstream myfile29;
	myfile29.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out29.str() + ".csv");
	myfile29.precision(25);

	t++;
	stringstream out30;
	out30 << t;
	ofstream myfile30;
	myfile30.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out30.str() + ".csv");
	myfile30.precision(25);

	t++;
	stringstream out31;
	out31 << t;
	ofstream myfile31;
	myfile31.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out31.str() + ".csv");
	myfile31.precision(25);

	t++;
	stringstream out32;
	out32 << t;
	ofstream myfile32;
	myfile32.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out32.str() + ".csv");
	myfile32.precision(25);

	t++;
	stringstream out33;
	out33 << t;
	ofstream myfile33;
	myfile33.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out33.str() + ".csv");
	myfile33.precision(25);

	t++;
	stringstream out34;
	out34 << t;
	ofstream myfile34;
	myfile34.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out34.str() + ".csv");
	myfile34.precision(25);

	t++;
	stringstream out35;
	out35 << t;
	ofstream myfile35;
	myfile35.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out35.str() + ".csv");
	myfile35.precision(25);

	t++;
	stringstream out36;
	out36 << t;
	ofstream myfile36;
	myfile36.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out36.str() + ".csv");
	myfile36.precision(25);

	t++;
	stringstream out37;
	out37 << t;
	ofstream myfile37;
	myfile37.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out37.str() + ".csv");
	myfile37.precision(25);

	t++;
	stringstream out38;
	out38 << t;
	ofstream myfile38;
	myfile38.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out38.str() + ".csv");
	myfile38.precision(25);

	t++;
	stringstream out39;
	out39 << t;
	ofstream myfile39;
	myfile39.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out39.str() + ".csv");
	myfile39.precision(25);

	t++;
	stringstream out40;
	out40 << t;
	ofstream myfile40;
	myfile40.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out40.str() + ".csv");
	myfile40.precision(25);

	t++;
	stringstream out41;
	out41 << t;
	ofstream myfile41;
	myfile41.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out41.str() + ".csv");
	myfile41.precision(25);

	t++;
	stringstream out42;
	out42 << t;
	ofstream myfile42;
	myfile42.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out42.str() + ".csv");
	myfile42.precision(25);

	t++;
	stringstream out43;
	out43 << t;
	ofstream myfile43;
	myfile43.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out43.str() + ".csv");
	myfile43.precision(25);

	t++;
	stringstream out44;
	out44 << t;
	ofstream myfile44;
	myfile44.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out44.str() + ".csv");
	myfile44.precision(25);

	t++;
	stringstream out45;
	out45 << t;
	ofstream myfile45;
	myfile45.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out45.str() + ".csv");
	myfile45.precision(25);

	t++;
	stringstream out46;
	out46 << t;
	ofstream myfile46;
	myfile46.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out46.str() + ".csv");
	myfile46.precision(25);

	t++;
	stringstream out47;
	out47 << t;
	ofstream myfile47;
	myfile47.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out47.str() + ".csv");
	myfile47.precision(25);

	t++;
	stringstream out48;
	out48 << t;
	ofstream myfile48;
	myfile48.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out48.str() + ".csv");
	myfile48.precision(25);

	t++;
	stringstream out49;
	out49 << t;
	ofstream myfile49;
	myfile49.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out49.str() + ".csv");
	myfile49.precision(25);

	t++;
	stringstream out50;
	out50 << t;
	ofstream myfile50;
	myfile50.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out50.str() + ".csv");
	myfile50.precision(25);

	t++;
	stringstream out51;
	out51 << t;
	ofstream myfile51;
	myfile51.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out51.str() + ".csv");
	myfile51.precision(25);

	t++;
	stringstream out52;
	out52 << t;
	ofstream myfile52;
	myfile52.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out52.str() + ".csv");
	myfile52.precision(25);

	t++;
	stringstream out53;
	out53 << t;
	ofstream myfile53;
	myfile53.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out53.str() + ".csv");
	myfile53.precision(25);

	t++;
	stringstream out54;
	out54 << t;
	ofstream myfile54;
	myfile54.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out54.str() + ".csv");
	myfile54.precision(25);

	t++;
	stringstream out55;
	out55 << t;
	ofstream myfile55;
	myfile55.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out55.str() + ".csv");
	myfile55.precision(25);

	t++;
	stringstream out56;
	out56 << t;
	ofstream myfile56;
	myfile56.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out56.str() + ".csv");
	myfile56.precision(25);

	t++;
	stringstream out57;
	out57 << t;
	ofstream myfile57;
	myfile57.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out57.str() + ".csv");
	myfile57.precision(25);

	t++;
	stringstream out58;
	out58 << t;
	ofstream myfile58;
	myfile58.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out58.str() + ".csv");
	myfile58.precision(25);

	t++;
	stringstream out59;
	out59 << t;
	ofstream myfile59;
	myfile59.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out59.str() + ".csv");
	myfile59.precision(25);

	t++;
	stringstream out60;
	out60 << t;
	ofstream myfile60;
	myfile60.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out60.str() + ".csv");
	myfile60.precision(25);

	t++;
	stringstream out61;
	out61 << t;
	ofstream myfile61;
	myfile61.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out61.str() + ".csv");
	myfile61.precision(25);

	t++;
	stringstream out62;
	out62 << t;
	ofstream myfile62;
	myfile62.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out62.str() + ".csv");
	myfile62.precision(25);

	t++;
	stringstream out63;
	out63 << t;
	ofstream myfile63;
	myfile63.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out63.str() + ".csv");
	myfile63.precision(25);

	t++;
	stringstream out64;
	out64 << t;
	ofstream myfile64;
	myfile64.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out64.str() + ".csv");
	myfile64.precision(25);

	t++;
	stringstream out65;
	out65 << t;
	ofstream myfile65;
	myfile65.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out65.str() + ".csv");
	myfile65.precision(25);

	t++;
	stringstream out66;
	out66 << t;
	ofstream myfile66;
	myfile66.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out66.str() + ".csv");
	myfile66.precision(25);

	t++;
	stringstream out67;
	out67 << t;
	ofstream myfile67;
	myfile67.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out67.str() + ".csv");
	myfile67.precision(25);

	t++;
	stringstream out68;
	out68 << t;
	ofstream myfile68;
	myfile68.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out68.str() + ".csv");
	myfile68.precision(25);

	t++;
	stringstream out69;
	out69 << t;
	ofstream myfile69;
	myfile69.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out69.str() + ".csv");
	myfile69.precision(25);

	t++;
	stringstream out70;
	out70 << t;
	ofstream myfile70;
	myfile70.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out70.str() + ".csv");
	myfile70.precision(25);

	t++;
	stringstream out71;
	out71 << t;
	ofstream myfile71;
	myfile71.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out71.str() + ".csv");
	myfile71.precision(25);

	t++;
	stringstream out72;
	out72 << t;
	ofstream myfile72;
	myfile72.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out72.str() + ".csv");
	myfile72.precision(25);

	t++;
	stringstream out73;
	out73 << t;
	ofstream myfile73;
	myfile73.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out73.str() + ".csv");
	myfile73.precision(25);

	t++;
	stringstream out74;
	out74 << t;
	ofstream myfile74;
	myfile74.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out74.str() + ".csv");
	myfile74.precision(25);

	t++;
	stringstream out75;
	out75 << t;
	ofstream myfile75;
	myfile75.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out75.str() + ".csv");
	myfile75.precision(25);

	t++;
	stringstream out76;
	out76 << t;
	ofstream myfile76;
	myfile76.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out76.str() + ".csv");
	myfile76.precision(25);

	t++;
	stringstream out77;
	out77 << t;
	ofstream myfile77;
	myfile77.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out77.str() + ".csv");
	myfile77.precision(25);

	t++;
	stringstream out78;
	out78 << t;
	ofstream myfile78;
	myfile78.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out78.str() + ".csv");
	myfile78.precision(25);

	t++;
	stringstream out79;
	out79 << t;
	ofstream myfile79;
	myfile79.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out79.str() + ".csv");
	myfile79.precision(25);

	t++;
	stringstream out80;
	out80 << t;
	ofstream myfile80;
	myfile80.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out80.str() + ".csv");
	myfile80.precision(25);

	t++;
	stringstream out81;
	out81 << t;
	ofstream myfile81;
	myfile81.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out81.str() + ".csv");
	myfile81.precision(25);

	t++;
	stringstream out82;
	out82 << t;
	ofstream myfile82;
	myfile82.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out82.str() + ".csv");
	myfile82.precision(25);

	t++;
	stringstream out83;
	out83 << t;
	ofstream myfile83;
	myfile83.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out83.str() + ".csv");
	myfile83.precision(25);

	t++;
	stringstream out84;
	out84 << t;
	ofstream myfile84;
	myfile84.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out84.str() + ".csv");
	myfile84.precision(25);

	t++;
	stringstream out85;
	out85 << t;
	ofstream myfile85;
	myfile85.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out85.str() + ".csv");
	myfile85.precision(25);

	t++;
	stringstream out86;
	out86 << t;
	ofstream myfile86;
	myfile86.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out86.str() + ".csv");
	myfile86.precision(25);

	t++;
	stringstream out87;
	out87 << t;
	ofstream myfile87;
	myfile87.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out87.str() + ".csv");
	myfile87.precision(25);

	t++;
	stringstream out88;
	out88 << t;
	ofstream myfile88;
	myfile88.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out88.str() + ".csv");
	myfile88.precision(25);

	t++;
	stringstream out89;
	out89 << t;
	ofstream myfile89;
	myfile89.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out89.str() + ".csv");
	myfile89.precision(25);

	t++;
	stringstream out90;
	out90 << t;
	ofstream myfile90;
	myfile90.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out90.str() + ".csv");
	myfile90.precision(25);

	t++;
	stringstream out91;
	out91 << t;
	ofstream myfile91;
	myfile91.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out91.str() + ".csv");
	myfile91.precision(25);

	t++;
	stringstream out92;
	out92 << t;
	ofstream myfile92;
	myfile92.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out92.str() + ".csv");
	myfile92.precision(25);

	t++;
	stringstream out93;
	out93 << t;
	ofstream myfile93;
	myfile93.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out93.str() + ".csv");
	myfile93.precision(25);

	t++;
	stringstream out94;
	out94 << t;
	ofstream myfile94;
	myfile94.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out94.str() + ".csv");
	myfile94.precision(25);

	t++;
	stringstream out95;
	out95 << t;
	ofstream myfile95;
	myfile95.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out95.str() + ".csv");
	myfile95.precision(25);

	t++;
	stringstream out96;
	out96 << t;
	ofstream myfile96;
	myfile96.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out96.str() + ".csv");
	myfile96.precision(25);

	t++;
	stringstream out97;
	out97 << t;
	ofstream myfile97;
	myfile97.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out97.str() + ".csv");
	myfile97.precision(25);

	t++;
	stringstream out98;
	out98 << t;
	ofstream myfile98;
	myfile98.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out98.str() + ".csv");
	myfile98.precision(25);

	t++;
	stringstream out99;
	out99 << t;
	ofstream myfile99;
	myfile99.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out99.str() + ".csv");
	myfile99.precision(25);

	t++;
	stringstream out100;
	out100 << t;
	ofstream myfile100;
	myfile100.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out100.str() + ".csv");
	myfile100.precision(25);

	t++;
	stringstream out101;
	out101 << t;
	ofstream myfile101;
	myfile101.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out101.str() + ".csv");
	myfile101.precision(25);

	t++;
	stringstream out102;
	out102 << t;
	ofstream myfile102;
	myfile102.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out102.str() + ".csv");
	myfile102.precision(25);

	t++;
	stringstream out103;
	out103 << t;
	ofstream myfile103;
	myfile103.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out103.str() + ".csv");
	myfile103.precision(25);

	t++;
	stringstream out104;
	out104 << t;
	ofstream myfile104;
	myfile104.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out104.str() + ".csv");
	myfile104.precision(25);

	t++;
	stringstream out105;
	out105 << t;
	ofstream myfile105;
	myfile105.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out105.str() + ".csv");
	myfile105.precision(25);

	t++;
	stringstream out106;
	out106 << t;
	ofstream myfile106;
	myfile106.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out106.str() + ".csv");
	myfile106.precision(25);

	t++;
	stringstream out107;
	out107 << t;
	ofstream myfile107;
	myfile107.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out107.str() + ".csv");
	myfile107.precision(25);

	t++;
	stringstream out108;
	out108 << t;
	ofstream myfile108;
	myfile108.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out108.str() + ".csv");
	myfile108.precision(25);

	t++;
	stringstream out109;
	out109 << t;
	ofstream myfile109;
	myfile109.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out109.str() + ".csv");
	myfile109.precision(25);

	t++;
	stringstream out110;
	out110 << t;
	ofstream myfile110;
	myfile110.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out110.str() + ".csv");
	myfile110.precision(25);

	t++;
	stringstream out111;
	out111 << t;
	ofstream myfile111;
	myfile111.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out111.str() + ".csv");
	myfile111.precision(25);

	t++;
	stringstream out112;
	out112 << t;
	ofstream myfile112;
	myfile112.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out112.str() + ".csv");
	myfile112.precision(25);

	t++;
	stringstream out113;
	out113 << t;
	ofstream myfile113;
	myfile113.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out113.str() + ".csv");
	myfile113.precision(25);

	t++;
	stringstream out114;
	out114 << t;
	ofstream myfile114;
	myfile114.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out114.str() + ".csv");
	myfile114.precision(25);

	t++;
	stringstream out115;
	out115 << t;
	ofstream myfile115;
	myfile115.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out115.str() + ".csv");
	myfile115.precision(25);

	t++;
	stringstream out116;
	out116 << t;
	ofstream myfile116;
	myfile116.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out116.str() + ".csv");
	myfile116.precision(25);

	t++;
	stringstream out117;
	out117 << t;
	ofstream myfile117;
	myfile117.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out117.str() + ".csv");
	myfile117.precision(25);

	t++;
	stringstream out118;
	out118 << t;
	ofstream myfile118;
	myfile118.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out118.str() + ".csv");
	myfile118.precision(25);

	t++;
	stringstream out119;
	out119 << t;
	ofstream myfile119;
	myfile119.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out119.str() + ".csv");
	myfile119.precision(25);

	t++;
	stringstream out120;
	out120 << t;
	ofstream myfile120;
	myfile120.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out120.str() + ".csv");
	myfile120.precision(25);

	t++;
	stringstream out121;
	out121 << t;
	ofstream myfile121;
	myfile121.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out121.str() + ".csv");
	myfile121.precision(25);

	t++;
	stringstream out122;
	out122 << t;
	ofstream myfile122;
	myfile122.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out122.str() + ".csv");
	myfile122.precision(25);

	t++;
	stringstream out123;
	out123 << t;
	ofstream myfile123;
	myfile123.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out123.str() + ".csv");
	myfile123.precision(25);

	t++;
	stringstream out124;
	out124 << t;
	ofstream myfile124;
	myfile124.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out124.str() + ".csv");
	myfile124.precision(25);

	t++;
	stringstream out125;
	out125 << t;
	ofstream myfile125;
	myfile125.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out125.str() + ".csv");
	myfile125.precision(25);

	t++;
	stringstream out126;
	out126 << t;
	ofstream myfile126;
	myfile126.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out126.str() + ".csv");
	myfile126.precision(25);

	t++;
	stringstream out127;
	out127 << t;
	ofstream myfile127;
	myfile127.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out127.str() + ".csv");
	myfile127.precision(25);

	t++;
	stringstream out128;
	out128 << t;
	ofstream myfile128;
	myfile128.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out128.str() + ".csv");
	myfile128.precision(25);

	t++;
	stringstream out129;
	out129 << t;
	ofstream myfile129;
	myfile129.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out129.str() + ".csv");
	myfile129.precision(25);

	t++;
	stringstream out130;
	out130 << t;
	ofstream myfile130;
	myfile130.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out130.str() + ".csv");
	myfile130.precision(25);

	t++;
	stringstream out131;
	out131 << t;
	ofstream myfile131;
	myfile131.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out131.str() + ".csv");
	myfile131.precision(25);

	t++;
	stringstream out132;
	out132 << t;
	ofstream myfile132;
	myfile132.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out132.str() + ".csv");
	myfile132.precision(25);

	t++;
	stringstream out133;
	out133 << t;
	ofstream myfile133;
	myfile133.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out133.str() + ".csv");
	myfile133.precision(25);

	t++;
	stringstream out134;
	out134 << t;
	ofstream myfile134;
	myfile134.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out134.str() + ".csv");
	myfile134.precision(25);

	t++;
	stringstream out135;
	out135 << t;
	ofstream myfile135;
	myfile135.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out135.str() + ".csv");
	myfile135.precision(25);

	t++;
	stringstream out136;
	out136 << t;
	ofstream myfile136;
	myfile136.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out136.str() + ".csv");
	myfile136.precision(25);

	t++;
	stringstream out137;
	out137 << t;
	ofstream myfile137;
	myfile137.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out137.str() + ".csv");
	myfile137.precision(25);

	t++;
	stringstream out138;
	out138 << t;
	ofstream myfile138;
	myfile138.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out138.str() + ".csv");
	myfile138.precision(25);

	t++;
	stringstream out139;
	out139 << t;
	ofstream myfile139;
	myfile139.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out139.str() + ".csv");
	myfile139.precision(25);

	t++;
	stringstream out140;
	out140 << t;
	ofstream myfile140;
	myfile140.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out140.str() + ".csv");
	myfile140.precision(25);

	t++;
	stringstream out141;
	out141 << t;
	ofstream myfile141;
	myfile141.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out141.str() + ".csv");
	myfile141.precision(25);

	t++;
	stringstream out142;
	out142 << t;
	ofstream myfile142;
	myfile142.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out142.str() + ".csv");
	myfile142.precision(25);

	t++;
	stringstream out143;
	out143 << t;
	ofstream myfile143;
	myfile143.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out143.str() + ".csv");
	myfile143.precision(25);

	t++;
	stringstream out144;
	out144 << t;
	ofstream myfile144;
	myfile144.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out144.str() + ".csv");
	myfile144.precision(25);

	t++;
	stringstream out145;
	out145 << t;
	ofstream myfile145;
	myfile145.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out145.str() + ".csv");
	myfile145.precision(25);

	t++;
	stringstream out146;
	out146 << t;
	ofstream myfile146;
	myfile146.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out146.str() + ".csv");
	myfile146.precision(25);

	t++;
	stringstream out147;
	out147 << t;
	ofstream myfile147;
	myfile147.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out147.str() + ".csv");
	myfile147.precision(25);

	t++;
	stringstream out148;
	out148 << t;
	ofstream myfile148;
	myfile148.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out148.str() + ".csv");
	myfile148.precision(25);

	t++;
	stringstream out149;
	out149 << t;
	ofstream myfile149;
	myfile149.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out149.str() + ".csv");
	myfile149.precision(25);

	t++;
	stringstream out150;
	out150 << t;
	ofstream myfile150;
	myfile150.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out150.str() + ".csv");
	myfile150.precision(25);

	t++;
	stringstream out151;
	out151 << t;
	ofstream myfile151;
	myfile151.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out151.str() + ".csv");
	myfile151.precision(25);

	t++;
	stringstream out152;
	out152 << t;
	ofstream myfile152;
	myfile152.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out152.str() + ".csv");
	myfile152.precision(25);

	t++;
	stringstream out153;
	out153 << t;
	ofstream myfile153;
	myfile153.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out153.str() + ".csv");
	myfile153.precision(25);

	t++;
	stringstream out154;
	out154 << t;
	ofstream myfile154;
	myfile154.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out154.str() + ".csv");
	myfile154.precision(25);

	t++;
	stringstream out155;
	out155 << t;
	ofstream myfile155;
	myfile155.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out155.str() + ".csv");
	myfile155.precision(25);

	t++;
	stringstream out156;
	out156 << t;
	ofstream myfile156;
	myfile156.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out156.str() + ".csv");
	myfile156.precision(25);

	t++;
	stringstream out157;
	out157 << t;
	ofstream myfile157;
	myfile157.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out157.str() + ".csv");
	myfile157.precision(25);

	t++;
	stringstream out158;
	out158 << t;
	ofstream myfile158;
	myfile158.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out158.str() + ".csv");
	myfile158.precision(25);

	t++;
	stringstream out159;
	out159 << t;
	ofstream myfile159;
	myfile159.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out159.str() + ".csv");
	myfile159.precision(25);

	t++;
	stringstream out160;
	out160 << t;
	ofstream myfile160;
	myfile160.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out160.str() + ".csv");
	myfile160.precision(25);

	t++;
	stringstream out161;
	out161 << t;
	ofstream myfile161;
	myfile161.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out161.str() + ".csv");
	myfile161.precision(25);

	t++;
	stringstream out162;
	out162 << t;
	ofstream myfile162;
	myfile162.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out162.str() + ".csv");
	myfile162.precision(25);

	t++;
	stringstream out163;
	out163 << t;
	ofstream myfile163;
	myfile163.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out163.str() + ".csv");
	myfile163.precision(25);

	t++;
	stringstream out164;
	out164 << t;
	ofstream myfile164;
	myfile164.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out164.str() + ".csv");
	myfile164.precision(25);

	t++;
	stringstream out165;
	out165 << t;
	ofstream myfile165;
	myfile165.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out165.str() + ".csv");
	myfile165.precision(25);

	t++;
	stringstream out166;
	out166 << t;
	ofstream myfile166;
	myfile166.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out166.str() + ".csv");
	myfile166.precision(25);

	t++;
	stringstream out167;
	out167 << t;
	ofstream myfile167;
	myfile167.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out167.str() + ".csv");
	myfile167.precision(25);

	t++;
	stringstream out168;
	out168 << t;
	ofstream myfile168;
	myfile168.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out168.str() + ".csv");
	myfile168.precision(25);

	t++;
	stringstream out169;
	out169 << t;
	ofstream myfile169;
	myfile169.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out169.str() + ".csv");
	myfile169.precision(25);

	t++;
	stringstream out170;
	out170 << t;
	ofstream myfile170;
	myfile170.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out170.str() + ".csv");
	myfile170.precision(25);

	t++;
	stringstream out171;
	out171 << t;
	ofstream myfile171;
	myfile171.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out171.str() + ".csv");
	myfile171.precision(25);

	t++;
	stringstream out172;
	out172 << t;
	ofstream myfile172;
	myfile172.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out172.str() + ".csv");
	myfile172.precision(25);

	t++;
	stringstream out173;
	out173 << t;
	ofstream myfile173;
	myfile173.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out173.str() + ".csv");
	myfile173.precision(25);

	t++;
	stringstream out174;
	out174 << t;
	ofstream myfile174;
	myfile174.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out174.str() + ".csv");
	myfile174.precision(25);

	t++;
	stringstream out175;
	out175 << t;
	ofstream myfile175;
	myfile175.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out175.str() + ".csv");
	myfile175.precision(25);

	t++;
	stringstream out176;
	out176 << t;
	ofstream myfile176;
	myfile176.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out176.str() + ".csv");
	myfile176.precision(25);

	t++;
	stringstream out177;
	out177 << t;
	ofstream myfile177;
	myfile177.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out177.str() + ".csv");
	myfile177.precision(25);

	t++;
	stringstream out178;
	out178 << t;
	ofstream myfile178;
	myfile178.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out178.str() + ".csv");
	myfile178.precision(25);

	t++;
	stringstream out179;
	out179 << t;
	ofstream myfile179;
	myfile179.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out179.str() + ".csv");
	myfile179.precision(25);

	t++;
	stringstream out180;
	out180 << t;
	ofstream myfile180;
	myfile180.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out180.str() + ".csv");
	myfile180.precision(25);

	t++;
	stringstream out181;
	out181 << t;
	ofstream myfile181;
	myfile181.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out181.str() + ".csv");
	myfile181.precision(25);

	t++;
	stringstream out182;
	out182 << t;
	ofstream myfile182;
	myfile182.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out182.str() + ".csv");
	myfile182.precision(25);

	t++;
	stringstream out183;
	out183 << t;
	ofstream myfile183;
	myfile183.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out183.str() + ".csv");
	myfile183.precision(25);

	t++;
	stringstream out184;
	out184 << t;
	ofstream myfile184;
	myfile184.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out184.str() + ".csv");
	myfile184.precision(25);

	t++;
	stringstream out185;
	out185 << t;
	ofstream myfile185;
	myfile185.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out185.str() + ".csv");
	myfile185.precision(25);

	t++;
	stringstream out186;
	out186 << t;
	ofstream myfile186;
	myfile186.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out186.str() + ".csv");
	myfile186.precision(25);

	t++;
	stringstream out187;
	out187 << t;
	ofstream myfile187;
	myfile187.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out187.str() + ".csv");
	myfile187.precision(25);

	t++;
	stringstream out188;
	out188 << t;
	ofstream myfile188;
	myfile188.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out188.str() + ".csv");
	myfile188.precision(25);

	t++;
	stringstream out189;
	out189 << t;
	ofstream myfile189;
	myfile189.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out189.str() + ".csv");
	myfile189.precision(25);

	t++;
	stringstream out190;
	out190 << t;
	ofstream myfile190;
	myfile190.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out190.str() + ".csv");
	myfile190.precision(25);

	t++;
	stringstream out191;
	out191 << t;
	ofstream myfile191;
	myfile191.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out191.str() + ".csv");
	myfile191.precision(25);

	t++;
	stringstream out192;
	out192 << t;
	ofstream myfile192;
	myfile192.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out192.str() + ".csv");
	myfile192.precision(25);

	t++;
	stringstream out193;
	out193 << t;
	ofstream myfile193;
	myfile193.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out193.str() + ".csv");
	myfile193.precision(25);

	t++;
	stringstream out194;
	out194 << t;
	ofstream myfile194;
	myfile194.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out194.str() + ".csv");
	myfile194.precision(25);

	t++;
	stringstream out195;
	out195 << t;
	ofstream myfile195;
	myfile195.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out195.str() + ".csv");
	myfile195.precision(25);

	t++;
	stringstream out196;
	out196 << t;
	ofstream myfile196;
	myfile196.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out196.str() + ".csv");
	myfile196.precision(25);

	t++;
	stringstream out197;
	out197 << t;
	ofstream myfile197;
	myfile197.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out197.str() + ".csv");
	myfile197.precision(25);

	t++;
	stringstream out198;
	out198 << t;
	ofstream myfile198;
	myfile198.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out198.str() + ".csv");
	myfile198.precision(25);

	t++;
	stringstream out199;
	out199 << t;
	ofstream myfile199;
	myfile199.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out199.str() + ".csv");
	myfile199.precision(25);

	t++;
	stringstream out200;
	out200 << t;
	ofstream myfile200;
	myfile200.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out200.str() + ".csv");
	myfile200.precision(25);

	t++;
	stringstream out201;
	out201 << t;
	ofstream myfile201;
	myfile201.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out201.str() + ".csv");
	myfile201.precision(25);

	t++;
	stringstream out202;
	out202 << t;
	ofstream myfile202;
	myfile202.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out202.str() + ".csv");
	myfile202.precision(25);

	t++;
	stringstream out203;
	out203 << t;
	ofstream myfile203;
	myfile203.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out203.str() + ".csv");
	myfile203.precision(25);

	t++;
	stringstream out204;
	out204 << t;
	ofstream myfile204;
	myfile204.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out204.str() + ".csv");
	myfile204.precision(25);

	t++;
	stringstream out205;
	out205 << t;
	ofstream myfile205;
	myfile205.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out205.str() + ".csv");
	myfile205.precision(25);

	/*new 205-231 */
	t++;
	stringstream out206;
	out206 << t;
	ofstream myfile206;
	myfile206.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out206.str() + ".csv");
	myfile206.precision(25);

	t++;
	stringstream out207;
	out207 << t;
	ofstream myfile207;
	myfile207.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out207.str() + ".csv");
	myfile207.precision(25);

	t++;
	stringstream out208;
	out208 << t;
	ofstream myfile208;
	myfile208.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out208.str() + ".csv");
	myfile208.precision(25);

	t++;
	stringstream out209;
	out209 << t;
	ofstream myfile209;
	myfile209.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out209.str() + ".csv");
	myfile209.precision(25);

	t++;
	stringstream out210;
	out210 << t;
	ofstream myfile210;
	myfile210.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out210.str() + ".csv");
	myfile210.precision(25);

	t++;
	stringstream out211;
	out211 << t;
	ofstream myfile211;
	myfile211.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out211.str() + ".csv");
	myfile211.precision(25);

	t++;
	stringstream out212;
	out212 << t;
	ofstream myfile212;
	myfile212.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out212.str() + ".csv");
	myfile212.precision(25);

	t++;
	stringstream out213;
	out213 << t;
	ofstream myfile213;
	myfile213.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out213.str() + ".csv");
	myfile213.precision(25);

	t++;
	stringstream out214;
	out214 << t;
	ofstream myfile214;
	myfile214.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out214.str() + ".csv");
	myfile214.precision(25);

	t++;
	stringstream out215;
	out215 << t;
	ofstream myfile215;
	myfile215.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out215.str() + ".csv");
	myfile215.precision(25);

	t++;
	stringstream out216;
	out216 << t;
	ofstream myfile216;
	myfile216.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out216.str() + ".csv");
	myfile216.precision(25);

	t++;
	stringstream out217;
	out217 << t;
	ofstream myfile217;
	myfile217.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out217.str() + ".csv");
	myfile217.precision(25);

	t++;
	stringstream out218;
	out218 << t;
	ofstream myfile218;
	myfile218.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out218.str() + ".csv");
	myfile218.precision(25);

	t++;
	stringstream out219;
	out219 << t;
	ofstream myfile219;
	myfile219.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out219.str() + ".csv");
	myfile219.precision(25);

	t++;
	stringstream out220;
	out220 << t;
	ofstream myfile220;
	myfile220.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out220.str() + ".csv");
	myfile220.precision(25);

	t++;
	stringstream out221;
	out221 << t;
	ofstream myfile221;
	myfile221.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out221.str() + ".csv");
	myfile221.precision(25);

	t++;
	stringstream out222;
	out222 << t;
	ofstream myfile222;
	myfile222.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out222.str() + ".csv");
	myfile222.precision(25);

	t++;
	stringstream out223;
	out223 << t;
	ofstream myfile223;
	myfile223.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out223.str() + ".csv");
	myfile223.precision(25);

	t++;
	stringstream out224;
	out224 << t;
	ofstream myfile224;
	myfile224.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out224.str() + ".csv");
	myfile224.precision(25);

	t++;
	stringstream out225;
	out225 << t;
	ofstream myfile225;
	myfile225.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out225.str() + ".csv");
	myfile225.precision(25);

	t++;
	stringstream out226;
	out226 << t;
	ofstream myfile226;
	myfile226.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out226.str() + ".csv");
	myfile226.precision(25);

	t++;
	stringstream out227;
	out227 << t;
	ofstream myfile227;
	myfile227.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out227.str() + ".csv");
	myfile227.precision(25);

	t++;
	stringstream out228;
	out228 << t;
	ofstream myfile228;
	myfile228.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out228.str() + ".csv");
	myfile228.precision(25);

	t++;
	stringstream out229;
	out229 << t;
	ofstream myfile229;
	myfile229.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out229.str() + ".csv");
	myfile229.precision(25);

	t++;
	stringstream out230;
	out230 << t;
	ofstream myfile230;
	myfile230.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out230.str() + ".csv");
	myfile230.precision(25);

	t++;
	stringstream out231;
	out231 << t;
	ofstream myfile231;
	myfile231.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\App_" + out231.str() + ".csv");
	myfile231.precision(25);

	
	string one;
	string two;
	string three;
	string four;
	string five;

	int user_count = 0;

	double total = 0;
	while (data.good())
	{

		total++;

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, '\n');
		
		

		/*  ---------------------------   */
		stringstream geek(one);
		double p1 = 0;
		geek >> p1;

		stringstream geek2(two);
		double p2 = 0;
		geek2 >> p2;

		stringstream geek3(three);
		double p3 = 0;
		geek3 >> p3;

		stringstream geek4(four);
		double p4 = 0;
		geek4 >> p4;

		stringstream geek5(five);
		double p5 = 0;
		geek5 >> p5;

		if (two =="yNf/R3X8fyXkOJm3ihXQcT0F52a8cDWPPRzTT6QFW8N+1QPfeKR5//6xyX0VYn7X")
		   myfile1 << p5 << ",";
		if (two == "gVb4X4iS13nJrM0KZsy7SrHzWAHix0CEPlK7/deV5vkwjt03xw5+QRtEJ9c8lHcD")
		   myfile2 << p5 << ",";
		if (two == "f8BtQHczrXFjnVxWl8Hqm1kH9UD/8nCvtXCLiFvVRvamqaxK0XRKUXz1ZOtcXQh6")
			myfile3 << p5 << ",";
		if (two == "8xjtVrJRJAyArNlbRgCftoNQOZiWl2eRP6uQavL9+6IULTcecLWn4PEn2c6kM9ha")
			myfile4 << p5 << ",";
		if (two == "zTE3f0H2n43tW+PA3OdonjUTWWxeyzF7xJk9QH9s/487J/EbT6DvZgVJ4AzkLBw5")
			myfile5 << p5 << ",";
		if (two == "oJXXVhFJaulRsMKo8iZ7PWewFKPbuwQtyXbR0IjQOWli8GZb66VGcdj+SpYd0Lsh")
			myfile6 << p5 << ",";
		if (two == "ef/Uph6sHGxOk9osPZNARutg0JeyAYagnPfsSr8TfYVf69sBLXydFO78zJLGpF/n")
			myfile7 << p5 << ",";
		if (two == "Dg/Eipen+htfY91WaoYpRA20XrBDDeUyUYjCeyJx/1T9Y0ge0MFQGch+aP/WiKIY")
			myfile8 << p5 << ",";
		if (two == "GEbO4GCYCUDDhA/p+yN1elnSN2p9oxmqwwdP/MiCKuZil6z/80BjHkiBn03rpspy")
			myfile9 << p5 << ",";
		if (two == "ZYtZQPMjjAwAO/1inl0hVnY1qbQTsVRzkw3sAAdDDN0avpZj8XVsmbZVJGdSMMJs")
			myfile10 << p5 << ",";
		if (two == "jU5bysNVd+WSutDiQwUTPwxhpFpNLe+4MVTLx19OXejsocU4UGdp062/jjmZHmJy")
			myfile11 << p5 << ",";
		if (two == "UB7m2LiuAwPW+9VBUgKarkZY2iU238Mct+KP+hxaPlXh7O5vIkd80sPnOuCMZJfH")
			myfile12 << p5 << ",";
		if (two == "U6ICBB2jjNegcsMyknapH8sHt5VOjtcFiryjPCZlJJBxk9Gl4UhW1A1PYpIwH2sr")
			myfile13 << p5 << ",";
		if (two == "yNf/R3X8fyXkOJm3ihXQcT0F52a8cDWPPRzTT6QFW8N+1QPfeKR5//6xyX0VYn7X")
			myfile14 << p5 << ",";
		if (two == "WulEtnpdd4Cb2vJKqNQQzufFjoxwEv1I/TsnJYS+RV2yOLeEN/dp9pqWymT1GlqP")
			myfile15 << p5 << ",";
		if (two == "I3WlW5a4rAwhh7gkdD4bFqdCpPlxuHJBovxCmLSQL/a8VjA3OdWso1s8IpVxiDrH")
			myfile16 << p5 << ",";
		if (two == "j9mDhLKahKDOkAAoeiXjv+JcJXAAKqJfRGtv/cxTha2FqZErbYfVgEd1oM6p+L6u")
			myfile17 << p5 << ",";
		if (two == "7pPCMOKxcej1YVIvrUGOYDzSliLLcSRdTnyF6Dtau6aIb+mR9zz+WBoaqvBP1EgM")
			myfile18 << p5 << ",";
		if (two == "zjBlanjORPr0Ku91HwIl+7Osa+6ApAkdUgrGSBqQDq5DyCQwK/ljfXKcuUMeWS/B")
			myfile19 << p5 << ",";
		if (two == "x4X9Xf0lJ81dq99ia7IP0+Fotz8leRbNj23585WIe/LX7Qa+Tjw0uF4ynk2oMTrU")
			myfile20 << p5 << ",";
		if (two == "Fs0RWRI+BTU8D0PDSTtkY/dKembMUBrUthPE/p8dkRXb+Eff7o88UGzmZ+v98Nx+")
			myfile21 << p5 << ",";
		if (two == "CN+SiZWbek6S/1iuQnGDuHxlnAO/ahD7u/e6bR0Ac4pMwdDZZhXMW1A60rSBTcWz")
			myfile22 << p5 << ",";
		if (two == "oY4PFnqfpD7zc6hXtAqone/fG9P96l7S3zGmxyFPhAxG3I5oQYWAsf5dBIt2iNKi")
			myfile23 << p5 << ",";
		if (two == "9KZrYg4+aCoZmHE0IzAYMbIMOmLTOJ1t+zL5BjLbRHl2FEhH8O3VEtbmnsdMcyGh")
			myfile24 << p5 << ",";
		if (two == "eN3W9pZnenLHR/1st/yBj8EYkfJc6fbSU5r/nbiEKuBTGXQWvLaXnhn7Ro5vqLyY")
			myfile25 << p5 << ",";
		if (two == "5KapxhMTJ0S+URdaKN/DGkRsF21iwxXyGrHB7f3CB20+nOzKDK8kQIYlLdWRYlKR")
			myfile26 << p5 << ",";
		if (two == "FXLAfeY7Xqos8x4+2DlVmcLdXrJZZ62JVDY3Fpkzen2lvJMxOvYtKVZlfVKmGJgR")
			myfile27 << p5 << ",";
		if (two == "mVa4KZ6O//zqNBfiJQ9/vStlPW4dQSPjiKwlrHSJYMpy92/7ENCqPyus3SMYUYZJ")
			myfile28 << p5 << ",";
		if (two == "rcqSFa1oo2BtGvSX9oblvbKUI6WmdLuThrZfMznkwNk4G17lzbStpQkumaY7rnrO")
			myfile29 << p5 << ",";
		if (two == "4Zs8R7YJjKZ6GPTUAWqp7rDRBY8weIeEafflEePF9ur/X+BE/4HfW6bAzVs5UfNr")
			myfile30 << p5 << ",";
		if (two == "Vvbj0IJzrk87L+EQhckyU53oLo+2tW+uHovl89ovDKL4uS3kLnvVbxOnVZiFrHc0")
			myfile31 << p5 << ",";
		if (two == "OcRL84edKZJZ9+vVUDuHsQHBRwM5HoFVoeYZ2fpZOFhxrDgHb/ci9HydyrNFBKE+")
			myfile32 << p5 << ",";
		if (two == "4BaYQCLoKM+sR22eeg0FQNG1SFFrAQpr/dIKP6VA+baPN8GxIoO5iWsg2mDz4jLa")
			myfile33 << p5 << ",";
		if (two == "6aYFGnH5DrZOQGDzKo2g+kc5GBqS9EI3L1cTEtzZD/AJW3pUvZr0DOBXVn6kar9T")
			myfile34 << p5 << ",";
		if (two == "A5MKq7BYw+qB83tDZ0h3IGwyVU5AWgi1ZHH6viwK4ONCAMSsiKrmT+VdUnyVRdQi")
			myfile35 << p5 << ",";
		if (two == "ZtqfHtAqSUukmdW6Ixj0QC+5dN9tnL27t73NOQfcgpa/QSWdBtEzz0ByrseaFts5")
			myfile36 << p5 << ",";
		if (two == "ecExvynUrLShF9XJua+Ho28AIvgipNenBFr81t9atoCM7LvfqmyxqJNYSUb9qW8J")
			myfile37 << p5 << ",";
		if (two == "uRmlLYBv7/eJM7fup3X3IjMnoUCc0QMi/8Iic75pYW30rT7J0DMFYa5IJC+sfBej")
			myfile38 << p5 << ",";
		if (two == "/qYUcpXfQMAxMRtJI3MHcg9+cUAheOSkXPLvIHcHi8bWWRfp0IAtmKe8SdLu0eT+")
			myfile39 << p5 << ",";
		if (two == "Bk6Krbp9OijvHJQyyltn5Y7TjRsXWzpwK7I2BkR0pQBjEkdXVSz0Acllwnzv/Ri6")
			myfile40 << p5 << ",";
		if (two == "hHo+r9tw45X27WSnBDjU/T9u4XFm+bpAg1G+yJRjAtVosr+Fm1Ynl+z9V5QvMqAj")
			myfile41 << p5 << ",";
		if (two == "up4yOvhjz+QMiX+vACyYmZgzRhIASBYvc8WcmNGHoLMszWnc1PH0qy4mzFe2wez4")
			myfile42 << p5 << ",";
		if (two == "bA2tJJsSFg7ZCsjBt1GRLiDzdqrM8po/vCebEme4QPhRKxb69KMCVL1Lljpu/3+9")
			myfile43 << p5 << ",";
		if (two == "lvOBEcngc/3Yt305BiFwrML327xFmSQi+QFwQZBAS5lgvKsqyLpMsYT0cZCaub6l")
			myfile44 << p5 << ",";
		if (two == "OvW1fP1VkpujEjZMBaVCd5Zj7DU2X7uS70QnTbaunpWxtIUzGLny1CwW2RrsZ/HP")
			myfile45 << p5 << ",";
		if (two == "X16z2cxWoWS2p8XrA5nFMSzqa+wGSnobXi7X9fG4QB+Qo9Btho1D4dgCJucg4D5l")
			myfile46 << p5 << ",";
		if (two == "pkIXjRx/GPQ7xU8jFREUciKuOOpE49JbzwnpK09gSCfOCzEIWT0oYHUYhI9jokkO")
			myfile47 << p5 << ",";
		if (two == "rYdTBu9KF5OzcKngM2ginb0Hxn+Cis8dCmkmvT1lo4cfkgIXIyDCB5kOtFW5yjDz")
			myfile48 << p5 << ",";
		if (two == "iHvKhF7Uk30Domr8uCORsSR904QSgslMQB+G+xe7JXNjm3QvGf2+VdR5gdk1mIrg")
			myfile49 << p5 << ",";
		if (two == "nEXGyIxF0UTklfB1oJsXDsUQ5os7rwjiul6bwqzoJJ6tbPjxNoWZXCDeF96Q6ynC")
			myfile50 << p5 << ",";

		/* New Apps*/
		if (two == "NL/R0ihfQrXNibeGXmHgGrBIDizDKGCmKMIJOHvLwXKjBVPEAJC73QDlRRJDe1KT")
			myfile51 << p5 << ",";
		if (two == "02ypkzq0vyvYVCqSQYjqIEkZgYZjddVdaZ2usBoKsK98k3KCHKPfRfYYuLYQj+I2")
			myfile52 << p5 << ",";
		if (two == "DV9CKTT9ElrfVj/Ow6Jkl1mUio4XBIHX0T9nzdAGU2ouVvuWiqkycR/HmAQ1aq+i")
			myfile53 << p5 << ",";
		if (two == "TGOkqLaoa8+9/gE3Qddi0VcV5fuMtz1DcDcZa9cvEuuQ94s/oJIWoQqylDfgKzCy")
			myfile54 << p5 << ",";
		if (two == "28P25Wf+gTpWS6mbRncA37MJcxE50+cDnKMG6PYVZvKuVw+M2lIY3u2wQ4qd7jam")
			myfile55 << p5 << ",";
		if (two == "cYugYImjey8hwJ3CMNfYAX/NyrxLfJlv3LE6apUpZEiknuAAL5VIXxIP92zqHSnJ")
			myfile56 << p5 << ",";
		if (two == "kcLNZyETdI1yFmaVl0c6xVJ3fifED1NYehwLvjea2mIqpgjpZb2K/PJTniuVjkGj")
			myfile57 << p5 << ",";
		if (two == "myqMf1mbYTgQYPziY0H36eNZfDtXUoe4VjZe2lDhN7eLhn5YXwq5E5DiR2uwykRu")
			myfile58 << p5 << ",";
		if (two == "faXNoLazRs508aNZmfWYj4ZlC/u0KFPWKoRR8mRjnvRnn+a+zrmlLvcPMzqlXt73")
			myfile59 << p5 << ",";
		if (two == "0SjR/C3z96Zb5gMGK2M87cK0hJ/5cslbK4rvLteCJmc2V21lsYE4JmkptFVJPlqw")
			myfile60 << p5 << ",";
		if (two == "oihJqA5/7N9jq5ewZFjj0D77Yb1MnG/shgS97Pmq/eEnwrUun50Xt8gVj5gMBoEA")
			myfile61 << p5 << ",";
		if (two == "Vk6Ir3hBspHY/ypz/UGauDFZ2jdm88wT5e8AkhsWFpxNFkV/FORfW7XC+sp8hhHd")
			myfile62 << p5 << ",";
		if (two == "C/+Lzeou2LizQySX67M2667gb5DYL8ydseiRleWccoX3IZo1wEXxWMFuLlHF23Qr")
			myfile63 << p5 << ",";
		if (two == "iZevtgAl0HI5hdm+uou75zIFDsD0BNXeEL68S8a0CaQA+yX4QMkd3LHvh185V5hs")
			myfile64 << p5 << ",";
		if (two == "LfQwIYEDJcs7vUACppNnvs5BQXT9F8GH74GHQIlJEkKHfyB+Kly20tdiYdHzlmbZ")
			myfile65 << p5 << ",";
		if (two == "WJ3tGObFv2PGf5pS+7+Hk6g83xwdE+SaPF6xEOo72E3iVjMtulv0BwKNe0aupaG6")
			myfile66 << p5 << ",";
		if (two == "TFS5vT3hLH4QhrkliyalzYhhBtw23fLIWffYgtQ58QSqC9nCBt7THzM4RL0GYKlF")
			myfile67 << p5 << ",";
		if (two == "5cABIv5rgzmu1wCUCVsd7w9NvSiKAOos+/T/HwUzMOWufpOQcVdkGv1Wj2YxCAF5")
			myfile68 << p5 << ",";
		if (two == "fZFMV7JihH/O0EFxE4KN+q0/ZE1r+7JH1WaUDQPHN0ZFG72v2Sx1KvAUJML+819z")
			myfile69 << p5 << ",";
		if (two == "VCyj8+WDGeawTbKHIoehGnPOQhB2b3UvGqj6jjR8LBplSDzGe1rkC26BQH9AbhW3")
			myfile70 << p5 << ",";
		if (two == "3PX7cSJGeh+nCHEImUfuoiXxv0ukvvQUQxRiqYK4qyuDlMru1tXKscbX1dOOp1CN")
			myfile71 << p5 << ",";
		if (two == "hvSoeHRu1PTHKwtSL0LxxiE4N/wM/xTf4+L61UL+bvWkxsERMu8MGgs6cbJFXyu/")
			myfile72 << p5 << ",";
		if (two == "ckZOxVGA15KlICneOKmROPL3iD7pLsQ9cOiOTPozo5ytHCQbXJ1PV54n4kTLE4Sj")
			myfile73 << p5 << ",";
		if (two == "chQrI6psarcrIXdXQZB4mHBPstwQBelmS3O2uFxAijkbsgGhEGF6LpBO4PndN5fd")
			myfile74 << p5 << ",";
		if (two == "JeaIN869z9QsEEgXO6ZYBSKCfvdTwSEb0qSCqTxFGH7a1pG6oTj6ca6m+XCxzJ2H")
			myfile75 << p5 << ",";
		if (two == "ZpMwnb9dWc2wJ8U8sM/xDN9n60fgii57F35vNtiUDPISjbJWDQ6wACH4th6A6GMu")
			myfile76 << p5 << ",";
		if (two == "TkxkQlPCxsGIQp1Cbt2Wh/tOGXc7KYM51TlJVnmIFAqOn10baeTFvVngFaiLAkDu")
			myfile77 << p5 << ",";
		if (two == "UXoXban2eQ7DpsWpkq0I/H21xRP6SqwlfGH/DxerBO1erClo97e+rhhnZP+dtKJ9")
			myfile78 << p5 << ",";
		if (two == "wBSAj80vcsficzeTavSm3kdA8ixLc2KzEK0d0tYfClyxmmEjiPe5BWkXMy7V1aOR")
			myfile79 << p5 << ",";
		if (two == "ZtH+o2LA+t2jsNxjSVmUXfEWBZ/K/QzmljtmaYpODfwDp9nZQbF6uQKgdm2mrKJ6")
			myfile80 << p5 << ",";
		if (two == "hxSdOXfdBNxacPfmxTjzoIwI0K15ddhbSdWSrFlQswhx/F6GbP4hTrHs5sadey0M")
			myfile81 << p5 << ",";
		if (two == "fDTzTl4YnPiwgiw65HwdBnwFW4BDdFkNkn/0Boh/DFrksZp+E/kzP0fnPWBlNJvx")
			myfile82 << p5 << ",";
		if (two == "J/J1cpRXim8NdkciZ5DzkulpoRr41Z0MNT8VuJT+ZNcyudhDdZjkITIBz2BFtVBf")
			myfile83 << p5 << ",";
		if (two == "zlfoPE6yz4qPBinPSpheS02PO79Q/xTrNd7L16P/pY3QAGV3znJQx9tpAC5FX3Cl")
			myfile84 << p5 << ",";
		if (two == "vIogsUicq/0LjmAJ4F4rR+ZQbI3bDpJzZnczbaS5Yo4fZX0tascK/4RxGmYbnYxc")
			myfile85 << p5 << ",";
		if (two == "gPLWtMtd1TKeLPjfserTlx1Fi+7XOmD1st2XhNS2pn/VO4rKCbKiaVcHEJLpq7pO")
			myfile86 << p5 << ",";
		if (two == "zss/AEEmzSQdCgh6NcNPRmPgbMbcrMfos3uU81fFY8FaGIbuPfPxZP+1XlxtBmzp")
			myfile87 << p5 << ",";
		if (two == "WN3C8R3g15gnD6zCZE9ZvWZTLTauUNwRYLE8fd0Tfg5WUv5xXYkil1X/luaZOPV2")
			myfile88 << p5 << ",";
		if (two == "VHZ5o2/T0s7RwYMFbL4YG5Io+eZVX5I1sDNuSMVTYUL/0Ohlv/IAhWjYm4sB4YSc")
			myfile89 << p5 << ",";
		if (two == "rBjH7WrJCeyxNVpZ5J9L6DfUYxWClx1wz6xMDrGoX2Hi1F2yHmjyOaTy0e3CTi8g")
			myfile90 << p5 << ",";
		if (two == "cguUNttemA0RxfBe5gnwQhG1cUi0DxRFMLnaa4/ISOsJNCa7RUFgTJSiDyS1Mbha")
			myfile91 << p5 << ",";
		if (two == "PpakhuuhrdIniBpgpUX8xqfuHlhdIGnKmvrv1GUQL7p5xSxgHMwmzrFP7eS1BBd9")
			myfile92 << p5 << ",";
		if (two == "njeoajdOAEjMFIEd3pKnHJ3naw+FnSmRAT19Tqg880f/dF6bTVnTYVkLLrLBmZYZ")
			myfile93 << p5 << ",";
		if (two == "r6CHLLsZ4/DiXAeZ7TojVWrGN5XDrsuj8bRLS8d8eQz1hyfx75p+n1UrYgmXrxEe")
			myfile94 << p5 << ",";
		if (two == "ZHCu/2qQsF/sWuSGR2itF88jKtwRy36nomPZqXt5827bCoWTItlq/5WGYCzbUuzN")
			myfile95 << p5 << ",";
		if (two == "CrmKrwlwGQ9En4rVFe5DXZpA3bL4PKCPFHwKWpovZKoqNo6UYdNgdcvMbfiZ5p0K")
			myfile96 << p5 << ",";
		if (two == "hbAbtRuLGTcndiJjOWLqoN5BhV1gvQ2guP4W/lRv9dGKHvivWHR+OkPiJFT9pVDq")
			myfile97 << p5 << ",";
		if (two == "7FT8KNiHgpdTjRs+baimn4HRlVBHHrF1WmJgUrgLsgJ9DSIzBc9UNlXTIWZzXY7i")
			myfile98 << p5 << ",";
		if (two == "ehseQuAs9vzo9tjJZVZWsjyIHITJ05XEaoF6KkDhGRxCoSs236CUbM6k6lKlQlCV")
			myfile99 << p5 << ",";
		if (two == "vHQ2NSLfJeSQOiQ9OPbqes71G0xhXJIBqcweOZfqing3TVPilNXk1ZM4X5dUaHuw")
			myfile100 << p5 << ",";
		/* New Apps*/
		if (two == "fTTTGF24gF4aKBg9oyZspi/PswZjAOdFR3Ex4F4CXldmcDdY70lcY9mXEvFwnROK")
			myfile101 << p5 << ",";
		if (two == "Bs6kes31gma4KiFtPkR6A/YWnYXpmLTwj4tfIFgOAE2wkrDo9QBBePI6g7lq49uJ")
			myfile102 << p5 << ",";
		if (two == "leaeqanOw3Cy3uN6H7tgCIbxWQX5ddahS5dgMjGWucgipSdj8QOHKygXb81DH4EC")
			myfile103 << p5 << ",";
		if (two == "guXfVyj7IYzNDTHS8A9C1iXleJxZcKQm4YiY9YY1lEEx11Q7Z9sM917JyIm2OA7h")
			myfile104 << p5 << ",";
		if (two == "oW8ocyHkkCBupu85tsFlwv2eOLF4Ffkm8lhg3kE/yoMsrKQGW6JIgRufx/PvyrJ7")
			myfile105 << p5 << ",";
		if (two == "lBe8nUqI8lmIxLZfd0XlZ/1GPDKMw/jFU8/DtW0qiCOfIVLAp+VqmG9zgHtzkFvb")
			myfile106 << p5 << ",";
		if (two == "2cXu+dEQtjw3Nj3Mo11tlobz6VDVYcGzq0SAKbceXkg0T0+FpPG2GVi8jXYlVOAl")
			myfile107 << p5 << ",";
		if (two == "plg52gEGIsfBZ4VcHMlaD5FC5iIhfJaTWmUL/R8kTqcWbOlEsNMxY1IXreGONCeE")
			myfile108 << p5 << ",";
		if (two == "XLZn6LkoqKohZfMsfqSVcxSDZqNhIM8aG2ugimQkhX1nw3odD1ztxTTTe1m/ic5u")
			myfile109 << p5 << ",";
		if (two == "BgFwBxavXoeyIWqsmjgt8Q+0A+n7z9Z3S/SHovqbTs4qCyuMqzUcPkuuk1rDq1Wj")
			myfile110 << p5 << ",";
		if (two == "Unf0DEjEEPum0kVnKo4s4EtTgf5pQ9Y8FwszKvmLqVqN7eJ8N1pf97WdTUP7BsGg")
			myfile111 << p5 << ",";
		if (two == "H94kKtIzBjenBPii6S3byu/GAue3fW2vhDxombe9AOBfYrDBKFaseiiDeHCT1Zpj")
			myfile112 << p5 << ",";
		if (two == "2zOV//yMuDrVsRYWdBxpu8WkkL5JtZzbxnIZbOLakZw/UjS4jdJ3ANmFgJFt5xbz")
			myfile113 << p5 << ",";
		if (two == "LyYN7J3LYjV3TE8U9llqbNtjPtG6pImzmnM9ygcSHOlxwHK8mLWYBVeo3VizzeED")
			myfile114 << p5 << ",";
		if (two == "g62Ztrkck0oFAWAiDvXOZFZ6usbx7bC3DjXb3fA5ZJ080WseRIQ4jNkbFDmTDeEt")
			myfile115 << p5 << ",";
		if (two == "jNO2pR51cSlkr3bcjVOfQS2KGQwPtxSMrw7ws9FYBCAdUeMeLTFJDnfCERUxBZNZ")
			myfile116 << p5 << ",";
		if (two == "1WBxiPX7yuf4c9da6x4KUFx7oB9OYnaF9UjXPUlQ/D57V73ui4PJ5gUfcMSsOCY/")
			myfile117 << p5 << ",";
		if (two == "gZxYPIKZnLtqggzhTzpg6j7O7QeIPyqzape8nuuYeB/a6P7xYxTIQZDpGRiUeFh0")
			myfile118 << p5 << ",";
		if (two == "9tnkECzlJqNYjUAHuWruT3pILBEC3IGpJIW5mr01w5wbmmpLR8P6zGnR5f2TKce/")
			myfile119 << p5 << ",";
		if (two == "dm2/wmfmt8KbPayyNF08HMhwRz5/zMs0wRtrBIur487j4UuJUpE8eY6HJqPgiUfj")
			myfile120 << p5 << ",";
		if (two == "p4Wep+HLOxkYUshdVBicO8SjW/x9gYybRlhxdQ1UiOLFeTyJltM80df/jKuP1jpC")
			myfile121 << p5 << ",";
		if (two == "dP49eFwPreYYpEpTzo0wKb5fSLg1zI7B1g6yTmuD6VZZ+lkP5n8OlgkkrMxm/DvY")
			myfile122 << p5 << ",";
		if (two == "SVXwLBvV/pieSXmEC037u1upuZRWetxRGZkDmrWmFLWmTSBIbtMQydTyUUNPUZ1/")
			myfile123 << p5 << ",";
		if (two == "yJe5ftz/a0MB494/elddvZ74IeStYQsVHNzEp2jd4pNbuC7RqsBaqyXWgvCdZa+q")
			myfile124 << p5 << ",";
		if (two == "Kh4QYpM7ohPG0c/7d27cGLBAbl0KwOmXPkfdKwNyg5yyEDmTPrjUIZgBGGSbO+Wf")
			myfile125 << p5 << ",";
		if (two == "SutBPn3KTKXGNxoCzJnIxL0FktSPVj88eGGgSOIv91iVFdcrTnvDkA1/BAtByWFE")
			myfile126 << p5 << ",";
		if (two == "Y/i5jUJrOjcxcQO9oyVoy7UB7ERlIqY1T0WfUojfzGMkAFhQiqqrGBtuJYj9TejM")
			myfile127 << p5 << ",";
		if (two == "6lS47hrW7jtifBhgZaWv9ekde/HD8bHF3g40+4e1uJbjjjV8j6Sf9gMjQOIk7TOQ")
			myfile128 << p5 << ",";
		if (two == "D9ZGzi+K0TsSpQfpx3sLnpYNjGxnQjA5+ahpKUBLKAqof1Uz3dl9450Y0dT+O/mY")
			myfile129 << p5 << ",";
		if (two == "X5SIZTmy8gyX3YSueNm1XRieuwIv+5hiSJdP4Btp0vnxdDuhaCEp9m2Ss9K4+pJv")
			myfile130 << p5 << ",";
		if (two == "sYg9r885ERdtvLwk8+DYodK91QsAsasaIAyuawG4VdYnRdshXw4WUbBzQXaFIVZh")
			myfile131 << p5 << ",";
		if (two == "GMRj97Sb5bcwYE/qrwY+gl+6vw2OGisogKICuVvXV2kiW0bukHs7M0XpsnpUrm5r")
			myfile132 << p5 << ",";
		if (two == "/sRJxBADGkYrL3kGufyF9pKfzi6T5YF9rU5slnvB5ryUul1+xZC3YQ2rEe9DkoLn")
			myfile133 << p5 << ",";
		if (two == "i325dXXRlQSK/dcxytSBcd6Qy0KrMjYm9qF0aPhq8+ZxRrFx76GJGs3WODytSa+v")
			myfile134 << p5 << ",";
		if (two == "nKHdCYJ1I3YEdhcPRe9tKFUHV1Q7eMfZBrms/J/RA704UhTC2IVYE/7saHoTBfIL")
			myfile135 << p5 << ",";
		if (two == "YnZbg2xWPCPpDIxmWpjXA1SamZhZ2h7OdILH6zF7nJ+YaLtr8eX59fYiMlQGXNRB")
			myfile136 << p5 << ",";
		if (two == "fUVXk3PDohHJH8ON/Jaj4ZNkXZ9uE9/1TooMyd7QA3Zgsz1VWpdDXMZoUtwAMm9M")
			myfile137 << p5 << ",";
		if (two == "ELUx4UZKNKvCOoAmZsILsC9Z3C2SpmM54XTyWkRV+SkWZVqwEuJ1xqizd8CIO+a5")
			myfile138 << p5 << ",";
		if (two == "TNQ7ssR90gsLNAOmbQcIn/PPNxptmiBavVdvgTZDqSrXEPrTOAcYl0Gg6gCJnpGH")
			myfile139 << p5 << ",";
		if (two == "EvncqblXtstU0Ip1pY6UkSCgZoFKZOc2iV0wgKEvvDVn/wbLoR4LBTb8Pq3Zskpm")
			myfile140 << p5 << ",";
		if (two == "qAxi1eYuV9CmCudewBejePCt30F+BXr5pwIxbpdupagdrWqpxCuRbpEJ2HW7OZHi")
			myfile141 << p5 << ",";
		if (two == "ll9OrpoTdmjz9tM3YypTdcbaxy6NV3+5wfGVgHbqFe2MEAGXCp1OpIk2LTv01ZtP")
			myfile142 << p5 << ",";
		if (two == "f5GVrlO9tcrNVBtO9uvzehvLCKiGeBEZBM0JTsRhxAGG6gjVN0BlK4HPIv+jzpcq")
			myfile143 << p5 << ",";
		if (two == "JNzYWeU5UXRzTb99YXPEBi1kTiHoVZJQ6hsPOWlvbW2IE8/jyZnFGMNwMdPexO22")
			myfile144 << p5 << ",";
		if (two == "V1Lpt70pOj5Y/jJ1CgZkQykumXYJ7Sv1yb+vkAQptlb9UhGQEHyhBIIdrdyeXx55")
			myfile145 << p5 << ",";
		if (two == "vGqiiqaGekaEnbNUix0ah3AgAWo/jMeKqJlaOTuEQyQ+B6ezSH0EyNRxlzskFP++")
			myfile146 << p5 << ",";
		if (two == "IBtikKpd2DU4aX4Eb/IAxGhFKIQvVMupQgQ3moLxpxzgGQu48FFyBONgoJhT7FRJ")
			myfile147 << p5 << ",";
		if (two == "tkK2nEV/zcKhPyDDeH0l+cHWGYn99rVO85rjwaoucGCQlZ1s3Tu5hkET7dsv46oC")
			myfile148 << p5 << ",";
		if (two == "DWoFQAL8pXYbOAq5BseTdqavUwxXKwJzkQvYBa8UjmhuZMtPROnhlN8Si+hZXyTj")
			myfile149 << p5 << ",";
		if (two == "VGKN/P/Sydrh13i9k0Y1j4daMg6CP1B24Ik/bacrFb29GLWxSxQrKIpUs0CwTkhY")
			myfile150 << p5 << ",";
		/* New Apps*/
		if (two == "VGcvRZuwhsOWg4z+UngA5AMaxzNvWZPajrfRSUMdKbKuW/xs1YiZeJa1svGw6GvZ")
			myfile151 << p5 << ",";
		if (two == "nRXMYQz5XUunjkeNPJObTotj7vpBrgRNHMA7swqeOyjARXvNadbrH+51SqZ0ztKH")
			myfile152 << p5 << ",";
		if (two == "MGe9pMuQq0maaWa98V+iK1eymYTtr5Gl5QlDBjQu27oUTB9ua8+5b51jOwn+dtEz")
			myfile153 << p5 << ",";
		if (two == "/f5dksGvuFfZtVBnEOw66+ZWbVP1cbZl4AqH/C22RLo3yuV6vQf01M9m6askrqg1")
			myfile154 << p5 << ",";
		if (two == "X+ntLerhHMJNkkWez6UmP+xbhehq4M0lr5YN7v5oeO13UktNap0XrnowCE2BMj/I")
			myfile155 << p5 << ",";
		if (two == "YW5w9yjMmDmm8kocE9Yuezi9C9F/j+JeC3bfJ0MjbWN1b3kDIFJaFZGW/Kfrvkey")
			myfile156 << p5 << ",";
		if (two == "5MJqsYqD6MYZt+C+ukeg+B/i9Z5/epMnEAATzP+YSb0V415ZHVlb4yHQqgqHsd4V")
			myfile157 << p5 << ",";
		if (two == "qL+zSq9u1JBuDg9Kdhbm+ymsh9wnUd9eyMhjpedmFP5PVO/C6dSJ7w2mevVStH32")
			myfile158 << p5 << ",";
		if (two == "RLnXJNnrE3NyXg9nvq1Ci33ypi4hE9Uy/ZhQ2ftteBfUJU0mSUHqxhOdy146uO6V")
			myfile159 << p5 << ",";
		if (two == "vzdDmBJhe688xnwJEr90L5NYFk7jGGm2H7TEQOCMfxnVPlBLQKHNsIOo/X9A4hbF")
			myfile160 << p5 << ",";
		if (two == "TQehyNWzqIefBUAXhiH2LOjPQz1APItmi0r0MZWYqXi2cblU8us1CrQpO0P5ikMV")
			myfile161 << p5 << ",";
		if (two == "UI/A1EWs2e0M7qwWeALfTi5wQtV7FEedwlNETsPMdzxASlve/TW4vZy4fzjrwngj")
			myfile162 << p5 << ",";
		if (two == "4nvtBXSfURR1iec6YIi/vv+QITn1b8yeKVfAVeAaki37AcA8/c5t7CsWDhCZgJyj")
			myfile163 << p5 << ",";
		if (two == "a/BNZOnvkO9X62zq6eRl4EouzXDxizas5T6M89U6/nbN8BLVNdF3SUFUK9OttmOP")
			myfile164 << p5 << ",";
		if (two == "AvQzeu6+Qg4AIaYhQbGUrqNUCSwaGFSK5Uv4VBzApHOZYsF9hFadPD8+0QH/V/ES")
			myfile165 << p5 << ",";
		if (two == "S2nOIYpbTR8hsrV1XntSNuhmJlU76DJmEVGqIj6plA0uZzB5JxNtiYcUjMpRHTdu")
			myfile166 << p5 << ",";
		if (two == "ylf2wZxBGywEpnUYNXv5K4V9lEhjYYz0UVY4GgzviAgukHCOcqX66vO/uiov1Jvm")
			myfile167 << p5 << ",";
		if (two == "xY8WpQDBuLH3c3U0Zv8o+iLIXdsnnxE3mxf7Pp8JnpmGP8IBeLQOTah9sPwRpUHM")
			myfile168 << p5 << ",";
		if (two == "o+DRJ+zL2G7hONeF8hQ3yx+7+t39UX1tALhuioROQCoZrUFDIVz959R7qE2Vx+Pq")
			myfile169 << p5 << ",";
		if (two == "Qi8QJ0yurvtZvrR1pk0PflfWOBXKcMS8QfDXVjFY9p6UQfKCh+R4YDOZQB60Nixh")
			myfile170 << p5 << ",";
		if (two == "X+MxRKzDJGq/arCZHP7uVU5EHnHsKnxa0F94JOQ87kSO8PL6AZYwGoYdUOoXXMMk")
			myfile171 << p5 << ",";
		if (two == "ouZYxM6rq28trmmHLFcNidZd2yX93OUsWI2TMkbuvgPUsmq26cjO5YglA+vc9pE1")
			myfile172 << p5 << ",";
		if (two == "vbMUToMvs83KL2rv0XsSHN5+AvnXOK97tk/6781g0BRAEsEjGV5F5rrHiC+kLyyB")
			myfile173 << p5 << ",";
		if (two == "KmAc2Hn3YJS9q3g648NnA+yjBEsVh5h777ZvWzRP86uRWM/5AZRgsvoGXKwKOp70")
			myfile174 << p5 << ",";
		if (two == "NB7ZLMwfTJyqaABC4njH9VpIQZflFLAEStabVgsPq/fDp3mYTizV18nJWTnkGYXk")
			myfile175 << p5 << ",";
		if (two == "gtfLEPsXoSrYCS/SmZkalr4x+oiAKN1Gggqmo5gNlp4jK7RwmmBKtJb57EkfBupt")
			myfile176 << p5 << ",";
		if (two == "vEtUNfikWcN3gQx8yz8Zruabh4QhMROhz7wSgfFMUVJVUI20/q54c2YFLzTve2gO")
			myfile177 << p5 << ",";
		if (two == "6loNog0hWzMxZ2jSqtgTCbmuYGb8ofghstp7P7b2LdWZ8V8Mm2ePywMc54J8B8Do")
			myfile178 << p5 << ",";
		if (two == "hCA+zCC/+rIr9ueh62xlzHhxLaWLo429Izqc9HbECdfQCBkL42vVFqOeupXsi10I")
			myfile179 << p5 << ",";
		if (two == "VHbHwzOPuC/y4Uw3Tv/shfIGZq2OLToxAufV6RIGh4NBydukPQc9ytkEyE5E5pqq")
			myfile180 << p5 << ",";
		if (two == "KX9jJtjno4B+cvRB+BNppVwG073n00+fiLGkzBCqbckaT1ZSklWBCSWirTwvftrV")
			myfile181 << p5 << ",";
		if (two == "RxLApsf2sRV/cXoY77mmXotKiZXEhPRnXt1gLUi66avt+as0yD6e4UGOIBfo3TAl")
			myfile182 << p5 << ",";
		if (two == "Du2+5VZVQPEAhmeUUc/AvqXoPWpeG8YsCHyC8lLdM+B1gDgX5ONNm3h1jHslDNFV")
			myfile183 << p5 << ",";
		if (two == "ToK4RPF7HgqqcdpCWT5Ew3pqjVfb6MzwpqFMD9Vw3Ye5lG8lyiiSIG5nHMjoXUmp")
			myfile184 << p5 << ",";
		if (two == "BqZDHYD8irARPRk9iTS2pGvKwSShAola68ymiOEzTINx98zHeys2eTcXPDh47ZDR")
			myfile185 << p5 << ",";
		if (two == "wko8RVXA00Y+qu3zJMCo4fXRppLFh9PMb95NQx3otRE+k5gVv473HZGbXUTJjiHs")
			myfile186 << p5 << ",";
		if (two == "/DC/fhHwZb8mE4meZb400SBBKXf3fJa/9zUxBUG9ye+wHLwCg0//fHN+C2jCxxxc")
			myfile187 << p5 << ",";
		if (two == "n2mVfhk/ouez8me3n8ajy9k+l6Z6bqhDxRPpT5JtNgcv/JjRqd0u+rTfAEoxF9Bu")
			myfile188 << p5 << ",";
		if (two == "ICuy+dnktp23pGsKAAsnuhKtet/DikbW9vOD4lImXIt1wh2hna+x9grvMYsQZMfP")
			myfile189 << p5 << ",";
		if (two == "DXSx/O2es/kbMIEtIiNWEPriftBSLrm924bXWtJ2uQQew/cDq2aLD9cUM4GaGUCQ")
			myfile190 << p5 << ",";
		if (two == "ocZgPFtObXCjcaQQX1SxLrycoHSX2lCW15lxzAjGPdVqguEufFMxUY/VGfrnTRzG")
			myfile191 << p5 << ",";
		if (two == "lWevj3virILbPLU8kVxLXIZ4hVYcKz/vUQAHPhcwBinm/FLeigSnNJfuzQ3j8XFb")
			myfile192 << p5 << ",";
		if (two == "mCyY0lEQSnm8+7XvzZ8VFAo0nXYzGyqUiE4xAHEVl5rJsIk+MFsRgcd8jDpvuOqV")
			myfile193 << p5 << ",";
		if (two == "gvLGWcpJytecikJmCgPyMW1ro2Qt+pvDY7Fj3dxhXzmFdi9g/O7kpKLJmHC5IaDi")
			myfile194 << p5 << ",";
		if (two == "WJ63tgCBBBEicbGTYufBfPOoHmv/ilzgpBKjK+c9el9H7T/HVHywljpqpiIEFrkY")
			myfile195 << p5 << ",";
		if (two == "YNJ9o2nkh54B0ZTKys0ce7My5DXhvRT7jeoty+prMB0sKcg1qelPMyNS0Dt4oF7W")
			myfile196 << p5 << ",";
		if (two == "z87OAWJyxVeQJ9ZNpm3my6OM5dGDgQDXC7Ks2kYL/O/hcSebg9wM3vngLWlt+Asf")
			myfile197 << p5 << ",";
		if (two == "2ldfXNeMwPeasOBL4QQR6xzZh78g3Y0zGm5lb+2w1mHnxQBnKl741cAab3BqjJSe")
			myfile198 << p5 << ",";
		if (two == "36CpLFm1MGYCuuN1j7bBmVWPZi6/9XI4VmK6s9oDGwHw+Z53LRUqu8tVmJTyaMFW")
			myfile199 << p5 << ",";
		if (two == "ddGx+/VOZh7cxTr1MYm2nGotc11VT1EogtsC7XK7qvyKB3BMJ9luJylOsRvbfr0m")
			myfile200 << p5 << ",";
		/* New Apps*/
		if (two == "YhMZM4qsnQqwEp++sIyIpo01PXC3OqPBRtZUJZlXBIc0g+IRzKyseMZHR1RcTc0k")
			myfile201 << p5 << ",";
		if (two == "MPq3W5RSiuyFIa26Dn9DJd7X3X/W8584oQ6DQhLm6V2Ku5IyALRLgokSQ+gAuHq5")
			myfile202 << p5 << ",";
		if (two == "DjltKt7pML/y3IQMPMS1hYo+FxLJCZD4PSykC0jPXJWvr8VnQ9yd6zYU/eXlyOz4")
			myfile203 << p5 << ",";
		if (two == "5Tv+e6Ptgwtlqxn5twHiIEiFmx6EC3TmOOaUVUB+3XZQlKYKOSsg3J+7QsbYYd2U")
			myfile204 << p5 << ",";
		if (two == "tqZ1p6Duh/xxLuqwsIXHkin2qgMY/7iDOc3200XLDARHriS195Kch3t2C9RTje9W")
			myfile205 << p5 << ",";
		if (two == "LlWsN1b+Ls5jt3q0x0fCoh2WtML9kTzqHWvA+fS8/PazreCvGj0sfFUSAznGXJRJ")
			myfile206 << p5 << ",";
		if (two == "/07OS/j+Flz2okPG8Ntpy6pg1zEadR7ROv6WvxKY04Smd1rOqDBSY+bXFA1H/525")
			myfile207 << p5 << ",";
		if (two == "NlYS1owpwUFtpASG5WmOCexjXbSEZLUD+6acZxzXW9kmc0oaZHI76B28HVy7RGgj")
			myfile208 << p5 << ",";
		if (two == "aH2kYHrhCURZg7pBudOH66yuNq2JVp0WUF3HxicX49mDw/mrccz0G1uWFSyuaRUs")
			myfile209 << p5 << ",";
		if (two == "A1rDLD+8Pj+LxaklIm4dib2mxpJVyINugZTAyO9Aa+Ce/I4L4LIJWL73PyqZnfOz")
			myfile210 << p5 << ",";
		if (two == "Ii3kBJ+TH85ApPCWPi/Hffb3GtOnr9u9yoAss0SNStbsVe9hg74O50CC8Jj+JryB")
			myfile211 << p5 << ",";
		if (two == "Dqnxsfh9GJk4qWlQwgfwFDR5Lyf5ysQkMvTimUWxut5EPbL3LSKcco8w89nOkA7O")
			myfile212 << p5 << ",";
		if (two == "FqYpILR/+thD+JFPG7C6IIjgkWmcmTQzLbbgSUG1FljsLk9fhzrex+kkH1qAXOWJ")
			myfile213 << p5 << ",";
		if (two == "6uo1RsK5K0iI7Y+fmcTGEBZg7tNlM3DvVOHFDkKQYkyDWXNnWT8DxlHRAZienpXY")
			myfile214 << p5 << ",";
		if (two == "wDuUSvsfFcXntU2x46hUMuTa3WVhWyYVUyVdRkFZVjH46OnZ4gT1HHPJ03gwAVHh")
			myfile215 << p5 << ",";
		if (two == "xWxi9smegw3J9T1oNOKaz3zOOPwZuBTtgba3jPEluLC28nEuedbJKZQLxIR4A3st")
			myfile216 << p5 << ",";
		if (two == "F5ZyHxj1Fl6YrRsGW1EFyFpMJN2vTW1iwjpbB5i7TUnBuyKyFLs/kvBFjYn7sEnK")
			myfile217 << p5 << ",";
		if (two == "JyNSyOj1vl1zjrr2UD8JqHtZgipKdjSPK3PMeVr9D3yuR3jR24YsVdwJMrD1BXtt")
			myfile218 << p5 << ",";
		if (two == "PN30OEiinNKHmjGrTW2CrTROeC3Dwk30Qa4agAuJfSrOUnPHfiVXN2RinlVsfU6/")
			myfile219 << p5 << ",";
		if (two == "epc3LPrsTmDeL/Wc3X/cECJnUbt1g3DPf+ImpCrCx3Hk2rRc9dgQlWLWvAci3z93")
			myfile220 << p5 << ",";
		if (two == "rvVjr/G0Mbp/J5A6rodh4sZRi1LkbuBrBcVxXCRADGCwDrnrRDzW1NTCF7A+S/JU")
			myfile221 << p5 << ",";
		if (two == "3mcq6nhWH3PqySLvlP16TpPR3hRDSFgr+xm50TmSwlR7A17aVVGzdAJCbNuQbF3p")
			myfile222 << p5 << ",";
		if (two == "ww8FdT9qTvViuwWo4KD8Mql9B4q+7n/WZDtAvCRF8xWwfuUPc2RSw4mP7p0bbENG")
			myfile223 << p5 << ",";
		if (two == "Ay4xfZBrAmTu8z9X2b5wtozC44GoxS0n5FxP+NGDV5+Hhln+7otKx+0VvUYkKlLI")
			myfile224 << p5 << ",";
		if (two == "uyOyeDCX1XnNw61/W2p0r9QH5IQzyqF70pJdsF0ZHjqC5OWaH3bXNN9/4Q3oRfGv")
			myfile225 << p5 << ",";
		if (two == "Cl1rMYhW9SCW1JaASotMMzH8MeHNoA9Hp38uL5GqTST3kEfp7VCY7b7rCmEHa3fz")
			myfile226 << p5 << ",";
		if (two == "5EdJSnR6MMVxLdl2ltSCiGYWeFhh/h9LQM/SDvU3tmCq+vrQ4P83bEXu4rBSE7/K")
			myfile227 << p5 << ",";
		if (two == "zAUob2BMaQptIUnaURED8jn3la57Jin7arm9x1NBD7TpdXXlpsBPLOxElWgr0TzL")
			myfile228 << p5 << ",";
		if (two == "5C40kzOMH3DIuToggPgdcdwmUs5t2UP+P3JR4AtG6kqtyqqutsReBGD8oVLN4x0k")
			myfile229 << p5 << ",";
		if (two == "RtqtPuL635zF55UlDUe2t28N3S3376Etbyqr45EI9EYA73KsmVQNKWV091hls+gB")
			myfile230 << p5 << ",";
		if (two == "0sR+C4cwv9ZIh2BNs9d2CZHNTpDE/tx7WIbJtI0WWs1oxyLmX4CUXnqi1Mj89v+n")
			myfile231 << p5 << ",";
		
  }

}

void SpecP()
{
	for (int i = 0;i<user_total;i++)
	{
		temporary_test[i] = new Users[20];
	}

	for (int t = 1;t <= 19;t++)
	{
		for (int i = 0;i < user_total;i++)
		{
			for (int h = 0;h < 20 - t;h++)
			{

				/* for predicting p */
				int sum4 = 0;
				for (int counter = 1; counter <= t;counter++)
					sum4 += (counter)*users[i][counter].p;
				sum4 /= ((t)*(t + 1) / 2);
				temporary_test[i][h].p = sum4;

			}
		}

		
		/*---------------------------------------*/
		stringstream out6;
		out6 << t;
		ofstream myfile6;
		myfile6.open("C:\\Users\\Erfan\\Desktop\\Processing Speed\\results_Sample_Data_50_p_" + out6.str() + ".csv");
		myfile6.precision(15);
		for (int i = 0;i < user_total;i++)
		{
			for (int n = 1;n <= t;n++)
				myfile6 << users[i][n].p << ",";
			for (int h = 0;h<20 - t;h++)
				myfile6 << temporary_test[i][h].p << ",";

			myfile6 << endl;
		}


	}

}

void ReadAfterBig()
{

	/* new dataset --- beginning-----  */
	for (int i = 0;i<user_total;i++)
	{
		users[i] = new Users[200];
		R[i] = new int[200];
		count_time[i] = 200;
	}

	cout.precision(15);

	ifstream data("C:\\Users\\Erfan\\Desktop\\New Dataset\\Sample_Data_210+21_b.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;

	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty,'\n');


		/*  ---------------------------   */
		stringstream geek(one);
		double b1 = 0;
		geek >> b1;
		users[user_count][1].md = b1*0.05;

		stringstream geek2(two);
		double b2 = 0;
		geek2 >> b2;
		users[user_count][2].md = b2*0.05;

		stringstream geek3(three);
		double b3 = 0;
		geek3 >> b3;
		users[user_count][3].md = b3*0.05;

		stringstream geek4(four);
		double b4 = 0;
		geek4 >> b4;
		users[user_count][4].md = b4*0.05;

		stringstream geek5(five);
		double b5 = 0;
		geek5 >> b5;
		users[user_count][5].md = b5*0.05;

		stringstream geek6(six);
		double b6 = 0;
		geek6 >> b6;
		users[user_count][6].md = b6*0.05;

		stringstream geek7(seven);
		double b7 = 0;
		geek7 >> b7;
		users[user_count][7].md = b7*0.05;

		stringstream geek8(eight);
		double b8 = 0;
		geek8 >> b8;
		users[user_count][8].md = b8*0.05;

		stringstream geek9(nine);
		double b9 = 0;
		geek9 >> b9;
		users[user_count][9].md = b9*0.05;

		stringstream geek10(ten);
		double b10 = 0;
		geek10 >> b10;
		users[user_count][10].md = b10*0.05;

		stringstream geek11(eleven);
		double b11 = 0;
		geek11 >> b11;
		users[user_count][11].md = b11*0.05;

		stringstream geek12(twelve);
		double b12 = 0;
		geek12 >> b12;
		users[user_count][12].md = b12*0.05;

		stringstream geek13(thirteen);
		double b13 = 0;
		geek13 >> b13;
		users[user_count][13].md = b13*0.05;

		stringstream geek14(fourteen);
		double b14 = 0;
		geek14 >> b14;
		users[user_count][14].md = b14*0.05;

		stringstream geek15(fifteen);
		double b15 = 0;
		geek15 >> b15;
		users[user_count][15].md = b15*0.05;

		stringstream geek16(sixteen);
		double b16 = 0;
		geek16 >> b16;
		users[user_count][16].md = b16*0.05;

		stringstream geek17(seventeen);
		double b17 = 0;
		geek17 >> b17;
		users[user_count][17].md = b17*0.05;

		stringstream geek18(eighteen);
		double b18 = 0;
		geek18 >> b18;
		users[user_count][18].md = b18*0.05;

		stringstream geek19(nineteen);
		double b19 = 0;
		geek19 >> b19;
		users[user_count][19].md = b19*0.05;

		stringstream geek20(twenty);
		double b20 = 0;
		geek20 >> b20;
		users[user_count][20].md = b20*0.05;

		user_count++;
	}

	stringstream out7;
	ofstream myfile7;
	myfile7.open("C:\\Users\\Erfan\\Desktop\\New Dataset\\Sample_Data_210+21_md.csv");
	myfile7.precision(15);
	for (int i = 0;i < user_total;i++)
	{
		for (int n = 1;n <= 20;n++)
			myfile7 << users[i][n].md << ",";
		myfile7 << endl;
	}
}

void MobilityConversion()
{
	const double halfC = M_PI / 180;

	double lat = 40.100523;
	double lon = -83.13209;
	double h = 274.93;

	double LatRad = lat*halfC;
	double LongRad = lon*halfC;

	double a = 6378.1370*1000;
	double b = 6356.7523*1000;
	double N;
	double e = 1 - pow(b, 2) / pow(a, 2);
	N = a / sqrt(1.0 - (e*pow(sin(LatRad), 2)));
	double cosLatRad = cos(LatRad);
	double cosLongiRad = cos(LongRad);
	double sinLatRad = sin(LatRad);
	double sinLongiRad = sin(LongRad);

	cout.precision(15);
	double x = (N + h)*cosLatRad*cosLongiRad;
	double y = (N + h)*cosLatRad*sinLongiRad;
	double z = ((pow(b, 2) /pow(a, 2))*N + h)*sinLatRad;

	cout<<"x "<<x<<endl;
	cout<<"y "<<y<<endl;
	cout<<"z "<<z<<endl;
}

void ForCoordinates()
{

	/* new dataset --- beginning-----  */
	for (int i = 0;i<user_total;i++)
	{
		users[i] = new Users[200];
		R[i] = new int[200];
		count_time[i] = 200;
	}


	cout.precision(15);

	string first;
	string second;
	string third;
	string fourth;
	string fifth;
	string sixth;
	string seventh;
	string eighth;
	string ninth;
	string tenth;
	string eleventh;

	ifstream data1("C:\\Users\\Erfan\\Desktop\\New Dataset\\New York.txt");

	if (!data1.is_open()) cout << "ERROR: File Open" << '\n';

	

	int user_count1 = 0;
	int index1 = 1;

	while (data1.good() && user_count1<user_total/3)
	{

		getline(data1, first, '\t');
		getline(data1, second, '\t');
		getline(data1, third, '\n');
		


		stringstream geek(second);
		double x = 0;
		geek >> x;
		cout.precision(second.size());

		stringstream geek2(third);
		double y = 0;
		geek2 >> y;
		cout.precision(third.size());


		int x1 = x - 6204.2827245510579814435914 + 588048.934405569;
		int y1 = y - 4386.3800635175866773352027 + (-4859258.6404757);

		if (x1 <= min_x)
			min_x = x1;
		if (x1>max_x)
			max_x = x1;

		if (y1 <= min_y)
			min_y = y1;
		if (y1>max_y)
			max_y = y1;


		if (user_count1 >= (user_total-21) / 3)
		{
			users[210 + user_count1 - (user_total - 21) / 3][index1].x = x1;
			users[210 + user_count1 - (user_total - 21) / 3][index1].y = y1;
			
		}
		else
		{
		users[user_count1][index1].x = x1;
		users[user_count1][index1].y = y1;
		}
		index1++;

		
		if (index1 == 21)
		{
			index1 = 1;
			user_count1++;
		}

	}

	/* start */
	ifstream data2("C:\\Users\\Erfan\\Desktop\\New Dataset\\Orlando.txt");

	if (!data2.is_open()) cout << "ERROR: File Open" << '\n';



	int user_count2 = 0;
	int index2 = 1;

	while (data2.good() && user_count2<user_total / 3)
	{

		getline(data2, first, '\t');
		getline(data2, second, '\t');
		getline(data2, third, '\n');

		stringstream geek(second);
		double x = 0;
		geek >> x;
		cout.precision(second.size());

		stringstream geek2(third);
		double y = 0;
		geek2 >> y;
		cout.precision(third.size());


		int x1 = x - (-57.456060041764366985717061) + 584240.202717688;
		int y1 = y - 169.9327558168544385353016 + (-4850678.11310361);

		if (x1 <= min_x)
			min_x = x1;
		if (x1>max_x)
			max_x = x1;

		if (y1 <= min_y)
			min_y = y1;
		if (y1>max_y)
			max_y = y1;

		if (user_count2 >= (user_total - 21) / 3)
		{
			users[217 + user_count2 - (user_total - 21) / 3][index2].x = x1;
			users[217 + user_count2 - (user_total - 21) / 3][index2].y = y1;

		}
		else
		{
			users[user_count2+70][index2].x = x1;
			users[user_count2+70][index2].y = y1;
		}
		index2++;


		if (index2 == 21)
		{
			index2 = 1;
			user_count2++;

		}

	}

	/*end */


	ifstream data3("C:\\Users\\Erfan\\Desktop\\DACT Strict-Dataset.csv");
	if (!data3.is_open()) cout << "ERROR: File Open" << '\n';


	int count3 = 0;
	int user_count3 = 0;
	int index3 = 1;

	while (data3.good() && user_count3<(user_total)/3)
	{

		getline(data3, first, ',');
		getline(data3, second, ',');
		getline(data3, third, ',');
		getline(data3, fourth, ',');
		getline(data3, fifth, ',');
		getline(data3, sixth, ',');
		getline(data3, seventh, ',');
		getline(data3, eighth, ',');
		getline(data3, ninth, ',');
		getline(data3, tenth, ',');

		count3++;

		stringstream geek(eighth);
		double x = 0;
		geek >> x;
		cout.precision(eighth.size());


		stringstream geek2(ninth);
		double y = 0;
		geek2 >> y;
		cout.precision(ninth.size());

		/* conversion to cartesian coordinates */
		const double halfC = M_PI / 180;

		double lat = x;
		double lon = y;
		double h = 274.93;

		double LatRad = lat*halfC;
		double LongRad = lon*halfC;

		double a2 = 6378.1370 * 1000;
		double b2 = 6356.7523 * 1000;
		double N;
		double e = 1 - pow(b2, 2) / pow(a2, 2);
		N = a2 / sqrt(1.0 - (e*pow(sin(LatRad), 2)));
		double cosLatRad = cos(LatRad);
		double cosLongiRad = cos(LongRad);
		double sinLatRad = sin(LatRad);
		double sinLongiRad = sin(LongRad);

		int x1 = (N + h)*cosLatRad*cosLongiRad;
		int y1 = (N + h)*cosLatRad*sinLongiRad;
		double z = ((pow(b2, 2) / pow(a2, 2))*N + h)*sinLatRad;


		if (x1 <= min_x)
			min_x = x1;
		if (x1>max_x)
			max_x = x1;

		if (y1 <= min_y)
			min_y = y1;
		if (y1>max_y)
			max_y = y1;

		if (count3 % 30 == 1)
		{
		    if(user_count3<(user_total-21)/3)
			{
			 users[user_count3+140][index3].x = x1;
			 users[user_count3+140][index3].y = y1;
			}
			else 
			{
				users[224 + user_count3 - (user_total - 21) / 3][index3].x = x1;
				users[224 + user_count3 - (user_total - 21) / 3][index3].y = y1;
			
			}
			index3++;
		}
		if (index3 == 21)
		{
			index3 = 1;
			user_count3++;
			count3 = 0;
		}

	}

	stringstream out6;
	ofstream myfile6;
	myfile6.open("C:\\Users\\Erfan\\Desktop\\New Dataset\\Sample_Data_210+21_x.csv");
	myfile6.precision(15);
	for (int i = 0;i < user_total;i++)
	{
		for (int n = 1;n <= 20;n++)
			myfile6 << users[i][n].x << ",";
		myfile6 << endl;
	}

	stringstream out7;
	ofstream myfile7;
	myfile7.open("C:\\Users\\Erfan\\Desktop\\New Dataset\\Sample_Data_210+21_y.csv");
	myfile7.precision(15);
	for (int i = 0;i < user_total;i++)
	{
		for (int n = 1;n <= 20;n++)
			myfile7 << users[i][n].y << ",";
		myfile7 << endl;
	}

}
void ReadMemory()
{
	/* new dataset --- beginning-----  */
	for (int i = 0;i<user_total;i++)
	{
		users[i] = new Users[200];
		R[i] = new int[200];
		count_time[i] = 200;
	}


	cout.precision(15);

	string first;
	string second;
	string third;
	string fourth;
	string fifth;
	string sixth;
	string seventh;
	string eighth;
	string ninth;
	string tenth;
	string eleventh;
	


	int Min_ID = INT_MAX;
	int Max_ID = INT_MIN;


	ifstream data("C:\\Users\\Erfan\\Desktop\\New Dataset\\rnd\\2013-9\\3.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';

	int index = 1;
	int count=0;
	int user_count = 0;

	while (data.good() && user_count<user_total)
	{

		getline(data, first, ';');
		getline(data, second, ';');
		getline(data, third, ';');
		getline(data, fourth, ';');
		getline(data, fifth, ';');
		getline(data, sixth, ';');
		getline(data, seventh, ';');
		getline(data, eighth, ';');
		getline(data, ninth, ';');
		getline(data, tenth, ';');
		getline(data, eleventh, '\n');
		count++;
		if (count > 1)
		{
			stringstream geek(seventh);
			double x = 0;
			geek >> x;
			cout.precision(seventh.size());

			users[user_count][index].id = x/1000;

			index++;
			if (index == 21)
			{
				index = 1;
				user_count++;
			}

			if(x<Min_ID)
			  Min_ID = x;
			if(x>Max_ID)
			  Max_ID = x;


		}


	}
	stringstream out7;
	ofstream myfile7;
	myfile7.open("C:\\Users\\Erfan\\Desktop\\New Dataset\\Sample_Data_210+21_id.csv");
	myfile7.precision(15);
	for (int i = 0;i < user_total;i++)
	{
		for (int n = 1;n <= 20;n++)
			myfile7 << users[i][n].id << ",";
		myfile7 << endl;
	}

	cout << "min ID " << Min_ID << endl;
	cout << "max MD " << Max_ID << endl;

}
void ReadBandwidth()
{
	/* new dataset --- beginning-----  */
	for (int i = 0;i<user_total;i++)
	{
		users[i] = new Users[200];
		R[i] = new int[200];
		count_time[i] = 200;
	}


	cout.precision(15);

	string first;
	string second;
	string third;
	string fourth;
	string fifth;
	string sixth;

	ifstream data("C:\\Users\\Erfan\\Desktop\\New Dataset\\Bandwidth.txt");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';

	cout.precision(25);

	int user_count = 0;
	int index = 1;

	while (data.good() && user_count<77)
	{

		getline(data, first, ' ');
		getline(data, second, ' ');
		getline(data, third, ' ');
		getline(data, fourth, ' ');
		getline(data, fifth, ' ');
		getline(data, sixth, '\n');


		stringstream geek(fifth);
		double x = 0;
		geek >> x;

		stringstream geek2(sixth);
		double y = 0;
		geek2 >> y;

		double temp = y/1000;
		double temp2 = (x*8*20)/1000000;
		int temp3 = temp2/temp;
		if (temp3 > 200)
		{

			users[user_count][index].b = temp2 / temp;
			index++;
		}


		if (index == 21)
		{
			index = 1;
			user_count++;
		}

	}

	stringstream out7;
	ofstream myfile7;
	myfile7.open("C:\\Users\\Erfan\\Desktop\\New Dataset\\Sample_Data_210+21_b.csv");
	myfile7.precision(15);
	for (int i = 0;i < user_total;i++)
	{
		for (int n = 1;n <= 20;n++)
			myfile7 << users[i][n].b << ",";
		myfile7 << endl;
	}
}
void DatasetGeneration()
{
	
    //ForCoordinates();  /* generating x and y */
	//ReadBig(); /* generating P */
	//ReadAfterBig(); /* read P and cobverting to MHz and integer */
	//ReadMemory(); /* generating id */;
	//ReadBandwidth(); /* generating bandwidth */
	ReadAfterBig();

}
void Results_MatrixCompletion()
{
	/* new dataset --- beginning-----  */
	for (int i = 0;i<user_total;i++)
	{

		First_Predict[i] = new Users[200];
		Second_Predict[i] = new Users[200];
		Third_Predict[i] = new Users[200];
		Fourth_Predict[i] = new Users[200];
		Fifth_Predict[i] =  new Users[200];

	}

	cout.precision(15);

  int t=1;
    
  while(t<=19)
  { 
	stringstream out;
	out << t;
	ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\210\\Results\\b\\results_Sample_Data_210_b_"+out.str()+".csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;


	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty, '\n');


		/*  ---------------------------   */
		stringstream geek(one);
		double b1 = 0;
		geek >> b1;

		stringstream geek2(two);
		double b2 = 0;
		geek2 >> b2;

		stringstream geek3(three);
		double b3 = 0;
		geek3 >> b3;

		stringstream geek4(four);
		double b4 = 0;
		geek4 >> b4;

		stringstream geek5(five);
		double b5 = 0;
		geek5 >> b5;

		stringstream geek6(six);
		double b6 = 0;
		geek6 >> b6;

		stringstream geek7(seven);
		double b7 = 0;
		geek7 >> b7;

		stringstream geek8(eight);
		double b8 = 0;
		geek8 >> b8;

		stringstream geek9(nine);
		double b9 = 0;
		geek9 >> b9;

		stringstream geek10(ten);
		double b10 = 0;
		geek10 >> b10;

		stringstream geek11(eleven);
		double b11 = 0;
		geek11 >> b11;

		stringstream geek12(twelve);
		double b12 = 0;
		geek12 >> b12;

		stringstream geek13(thirteen);
		double b13 = 0;
		geek13 >> b13;

		stringstream geek14(fourteen);
		double b14 = 0;
		geek14 >> b14;

		stringstream geek15(fifteen);
		double b15 = 0;
		geek15 >> b15;

		stringstream geek16(sixteen);
		double b16 = 0;
		geek16 >> b16;

		stringstream geek17(seventeen);
		double b17 = 0;
		geek17 >> b17;

		stringstream geek18(eighteen);
		double b18 = 0;
		geek18 >> b18;

		stringstream geek19(nineteen);
		double b19 = 0;
		geek19 >> b19;

		stringstream geek20(twenty);
		double b20 = 0;
		geek20 >> b20;

		if (t == 1)
		{
		  First_Predict[user_count][t].b = b2;
		  Second_Predict[user_count][t].b = b3;
		  Third_Predict[user_count][t].b = b4;
		  Fourth_Predict[user_count][t].b = b5;
		  Fifth_Predict[user_count][t].b = b6;
		}
		if (t == 2)
		{
			First_Predict[user_count][t].b = b3;
			Second_Predict[user_count][t].b = b4;
			Third_Predict[user_count][t].b = b5;
			Fourth_Predict[user_count][t].b = b6;
			Fifth_Predict[user_count][t].b = b7;
		}
		if (t == 3)
		{
			First_Predict[user_count][t].b = b4;
			Second_Predict[user_count][t].b = b5;
			Third_Predict[user_count][t].b = b6;
			Fourth_Predict[user_count][t].b = b7;
			Fifth_Predict[user_count][t].b = b8;
		}
		if (t == 4)
		{
			First_Predict[user_count][t].b = b5;
			Second_Predict[user_count][t].b = b6;
			Third_Predict[user_count][t].b = b7;
			Fourth_Predict[user_count][t].b = b8;
			Fifth_Predict[user_count][t].b = b9;
		}
		if (t == 5)
		{
			First_Predict[user_count][t].b = b6;
			Second_Predict[user_count][t].b = b7;
			Third_Predict[user_count][t].b = b8;
			Fourth_Predict[user_count][t].b = b9;
			Fifth_Predict[user_count][t].b = b10;
		}
		if (t == 6)
		{
			First_Predict[user_count][t].b = b7;
			Second_Predict[user_count][t].b = b8;
			Third_Predict[user_count][t].b = b9;
			Fourth_Predict[user_count][t].b = b10;
			Fifth_Predict[user_count][t].b = b11;
		}
		if (t == 7)
		{
			First_Predict[user_count][t].b = b8;
			Second_Predict[user_count][t].b = b9;
			Third_Predict[user_count][t].b = b10;
			Fourth_Predict[user_count][t].b = b11;
			Fifth_Predict[user_count][t].b = b12;
		}
		if (t == 8)
		{
			First_Predict[user_count][t].b = b9;
			Second_Predict[user_count][t].b = b10;
			Third_Predict[user_count][t].b = b11;
			Fourth_Predict[user_count][t].b = b12;
			Fifth_Predict[user_count][t].b = b13;
		}
		if (t == 9)
		{
			First_Predict[user_count][t].b = b10;
			Second_Predict[user_count][t].b = b11;
			Third_Predict[user_count][t].b = b12;
			Fourth_Predict[user_count][t].b = b13;
			Fifth_Predict[user_count][t].b = b14;
		}
		if (t == 10)
		{
			First_Predict[user_count][t].b = b11;
			Second_Predict[user_count][t].b = b12;
			Third_Predict[user_count][t].b = b13;
			Fourth_Predict[user_count][t].b = b14;
			Fifth_Predict[user_count][t].b = b15;
		}
		if (t == 11)
		{
			First_Predict[user_count][t].b = b12;
			Second_Predict[user_count][t].b = b13;
			Third_Predict[user_count][t].b = b14;
			Fourth_Predict[user_count][t].b = b15;
			Fifth_Predict[user_count][t].b = b16;
		}
		if (t == 12)
		{
			First_Predict[user_count][t].b = b13;
			Second_Predict[user_count][t].b = b14;
			Third_Predict[user_count][t].b = b15;
			Fourth_Predict[user_count][t].b = b16;
			Fifth_Predict[user_count][t].b = b17;
		}
		if (t == 13)
		{
			First_Predict[user_count][t].b = b14;
			Second_Predict[user_count][t].b = b15;
			Third_Predict[user_count][t].b = b16;
			Fourth_Predict[user_count][t].b = b17;
			Fifth_Predict[user_count][t].b = b18;
		}
		if (t == 14)
		{
			First_Predict[user_count][t].b = b15;
			Second_Predict[user_count][t].b = b16;
			Third_Predict[user_count][t].b = b17;
			Fourth_Predict[user_count][t].b = b18;
			Fifth_Predict[user_count][t].b = b19;
		}
		if (t == 15)
		{
			First_Predict[user_count][t].b = b16;
			Second_Predict[user_count][t].b = b17;
			Third_Predict[user_count][t].b = b18;
			Fourth_Predict[user_count][t].b = b19;
			Fifth_Predict[user_count][t].b = b20;
		}
		if (t == 16)
		{
			First_Predict[user_count][t].b = b17;
			Second_Predict[user_count][t].b = b18;
			Third_Predict[user_count][t].b = b19;
			Fourth_Predict[user_count][t].b = b20;
		}
		if (t == 17)
		{
			First_Predict[user_count][t].b = b18;
			Second_Predict[user_count][t].b = b19;
			Third_Predict[user_count][t].b = b20;

		}
		if (t == 18)
		{
			First_Predict[user_count][t].b = b19;
			Second_Predict[user_count][t].b = b20;
		}
		if (t == 19)
		{
			First_Predict[user_count][t].b = b20;
		}
		if (t == 20)
		{
			
		}

		user_count++;
	}
	t++;
   }

   t=1;
   /* id */
   while (t <= 19)
   {
	   stringstream out;
	   out << t;
	   ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\210\\Results\\id\\results_Sample_Data_210_id_" + out.str() + ".csv");

	   if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	   string one;
	   string two;
	   string three;
	   string four;
	   string five;
	   string six;
	   string seven;
	   string eight;
	   string nine;
	   string ten;
	   string eleven;
	   string twelve;
	   string thirteen;
	   string fourteen;
	   string fifteen;
	   string sixteen;
	   string seventeen;
	   string eighteen;
	   string nineteen;
	   string twenty;

	   int user_count = 0;

	   while (data.good() && user_count<user_total)
	   {

		   getline(data, one, ',');
		   getline(data, two, ',');
		   getline(data, three, ',');
		   getline(data, four, ',');
		   getline(data, five, ',');
		   getline(data, six, ',');
		   getline(data, seven, ',');
		   getline(data, eight, ',');
		   getline(data, nine, ',');
		   getline(data, ten, ',');
		   getline(data, eleven, ',');
		   getline(data, twelve, ',');
		   getline(data, thirteen, ',');
		   getline(data, fourteen, ',');
		   getline(data, fifteen, ',');
		   getline(data, sixteen, ',');
		   getline(data, seventeen, ',');
		   getline(data, eighteen, ',');
		   getline(data, nineteen, ',');
		   getline(data, twenty, '\n');


		   /*  ---------------------------   */
		   stringstream geek(one);
		   double id1 = 0;
		   geek >> id1;

		   stringstream geek2(two);
		   double id2 = 0;
		   geek2 >> id2;

		   stringstream geek3(three);
		   double id3 = 0;
		   geek3 >> id3;

		   stringstream geek4(four);
		   double id4 = 0;
		   geek4 >> id4;

		   stringstream geek5(five);
		   double id5 = 0;
		   geek5 >> id5;

		   stringstream geek6(six);
		   double id6 = 0;
		   geek6 >> id6;

		   stringstream geek7(seven);
		   double id7 = 0;
		   geek7 >> id7;

		   stringstream geek8(eight);
		   double id8 = 0;
		   geek8 >> id8;

		   stringstream geek9(nine);
		   double id9 = 0;
		   geek9 >> id9;

		   stringstream geek10(ten);
		   double id10 = 0;
		   geek10 >> id10;

		   stringstream geek11(eleven);
		   double id11 = 0;
		   geek11 >> id11;

		   stringstream geek12(twelve);
		   double id12 = 0;
		   geek12 >> id12;

		   stringstream geek13(thirteen);
		   double id13 = 0;
		   geek13 >> id13;

		   stringstream geek14(fourteen);
		   double id14 = 0;
		   geek14 >> id14;

		   stringstream geek15(fifteen);
		   double id15 = 0;
		   geek15 >> id15;

		   stringstream geek16(sixteen);
		   double id16 = 0;
		   geek16 >> id16;

		   stringstream geek17(seventeen);
		   double id17 = 0;
		   geek17 >> id17;

		   stringstream geek18(eighteen);
		   double id18 = 0;
		   geek18 >> id18;

		   stringstream geek19(nineteen);
		   double id19 = 0;
		   geek19 >> id19;

		   stringstream geek20(twenty);
		   double id20 = 0;
		   geek20 >> id20;

		   if (t == 1)
		   {
			   First_Predict[user_count][t].id = id2*100;
			   Second_Predict[user_count][t].id = id3*100;
			   Third_Predict[user_count][t].id = id4*100;
			   Fourth_Predict[user_count][t].id = id5*100;
			   Fifth_Predict[user_count][t].id = id6*100;
		   }
		   if (t == 2)
		   {
			   First_Predict[user_count][t].id = id3*100;
			   Second_Predict[user_count][t].id = id4*100;
			   Third_Predict[user_count][t].id = id5*100;
			   Fourth_Predict[user_count][t].id = id6*100;
			   Fifth_Predict[user_count][t].id = id7*100;
		   }
		   if (t == 3)
		   {
			   First_Predict[user_count][t].id = id4*100;
			   Second_Predict[user_count][t].id = id5*100;
			   Third_Predict[user_count][t].id = id6*100;
			   Fourth_Predict[user_count][t].id = id7*100;
			   Fifth_Predict[user_count][t].id = id8*100;
		   }
		   if (t == 4)
		   {
			   First_Predict[user_count][t].id = id5*100;
			   Second_Predict[user_count][t].id = id6*100;
			   Third_Predict[user_count][t].id = id7*100;
			   Fourth_Predict[user_count][t].id = id8*100;
			   Fifth_Predict[user_count][t].id = id9*100;
		   }
		   if (t == 5)
		   {
			   First_Predict[user_count][t].id = id6*100;
			   Second_Predict[user_count][t].id = id7*100;
			   Third_Predict[user_count][t].id = id8*100;
			   Fourth_Predict[user_count][t].id = id9*100;
			   Fifth_Predict[user_count][t].id = id10*100;
		   }
		   if (t == 6)
		   {
			   First_Predict[user_count][t].id = id7*100;
			   Second_Predict[user_count][t].id = id8*100;
			   Third_Predict[user_count][t].id = id9*100;
			   Fourth_Predict[user_count][t].id = id10*100;
			   Fifth_Predict[user_count][t].id = id11*100;
		   }
		   if (t == 7)
		   {
			   First_Predict[user_count][t].id = id8*100;
			   Second_Predict[user_count][t].id = id9*100;
			   Third_Predict[user_count][t].id = id10*100;
			   Fourth_Predict[user_count][t].id = id11*100;
			   Fifth_Predict[user_count][t].id = id12*100;
		   }
		   if (t == 8)
		   {
			   First_Predict[user_count][t].id = id9*100;
			   Second_Predict[user_count][t].id = id10*100;
			   Third_Predict[user_count][t].id = id11*100;
			   Fourth_Predict[user_count][t].id = id12*100;
			   Fifth_Predict[user_count][t].id = id13*100;
		   }
		   if (t == 9)
		   {
			   First_Predict[user_count][t].id = id10*100;
			   Second_Predict[user_count][t].id = id11*100;
			   Third_Predict[user_count][t].id = id12*100;
			   Fourth_Predict[user_count][t].id = id13*100;
			   Fifth_Predict[user_count][t].id = id14*100;
		   }
		   if (t == 10)
		   {
			   First_Predict[user_count][t].id = id11*100;
			   Second_Predict[user_count][t].id = id12*100;
			   Third_Predict[user_count][t].id = id13*100;
			   Fourth_Predict[user_count][t].id = id14*100;
			   Fifth_Predict[user_count][t].id = id15*100;
		   }
		   if (t == 11)
		   {
			   First_Predict[user_count][t].id = id12*100;
			   Second_Predict[user_count][t].id = id13*100;
			   Third_Predict[user_count][t].id = id14*100;
			   Fourth_Predict[user_count][t].id = id15*100;
			   Fifth_Predict[user_count][t].id = id16*100;
		   }
		   if (t == 12)
		   {
			   First_Predict[user_count][t].id = id13*100;
			   Second_Predict[user_count][t].id = id14*100;
			   Third_Predict[user_count][t].id = id15*100;
			   Fourth_Predict[user_count][t].id = id16*100;
			   Fifth_Predict[user_count][t].id = id17*100;
		   }
		   if (t == 13)
		   {
			   First_Predict[user_count][t].id = id14*100;
			   Second_Predict[user_count][t].id = id15*100;
			   Third_Predict[user_count][t].id = id16*100;
			   Fourth_Predict[user_count][t].id = id17*100;
			   Fifth_Predict[user_count][t].id = id18*100;
		   }
		   if (t == 14)
		   {
			   First_Predict[user_count][t].id = id15*100;
			   Second_Predict[user_count][t].id = id16*100;
			   Third_Predict[user_count][t].id = id17*100;
			   Fourth_Predict[user_count][t].id = id18*100;
			   Fifth_Predict[user_count][t].id = id19*100;
		   }
		   if (t == 15)
		   {
			   First_Predict[user_count][t].id = id16*100;
			   Second_Predict[user_count][t].id = id17*100;
			   Third_Predict[user_count][t].id = id18*100;
			   Fourth_Predict[user_count][t].id = id19*100;
			   Fifth_Predict[user_count][t].id = id20*100;
		   }
		   if (t == 16)
		   {
			   First_Predict[user_count][t].id = id17*100;
			   Second_Predict[user_count][t].id = id18*100;
			   Third_Predict[user_count][t].id = id19*100;
			   Fourth_Predict[user_count][t].id = id20*100;
		   }
		   if (t == 17)
		   {
			   First_Predict[user_count][t].id = id18*100;
			   Second_Predict[user_count][t].id = id19*100;
			   Third_Predict[user_count][t].id = id20*100;
		   }
		   if (t == 18)
		   {
			   First_Predict[user_count][t].id = id19*100;
			   Second_Predict[user_count][t].id = id20*100;
		   }
		   if (t == 19)
		   {
			   First_Predict[user_count][t].id = id20*100;
		   }
		   if (t == 20)
		   {
			   //First_Predict[user_count][t].b = b2;
			   //Second_Predict[user_count][t].b = b3;
			   //Third_Predict[user_count][t].b = b4;
			   //Fourth_Predict[user_count][t].b = b5;
			   //Fifth_Predict[user_count][t].b = b6;
		   }

		   user_count++;
	   }
	   t++;
   }

   /* md */
   t = 1;

   while (t <= 19)
   {
	   stringstream out;
	   out << t;
	   ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\210\\Results\\md\\results_Sample_Data_210_md_" + out.str() + ".csv");

	   if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	   string one;
	   string two;
	   string three;
	   string four;
	   string five;
	   string six;
	   string seven;
	   string eight;
	   string nine;
	   string ten;
	   string eleven;
	   string twelve;
	   string thirteen;
	   string fourteen;
	   string fifteen;
	   string sixteen;
	   string seventeen;
	   string eighteen;
	   string nineteen;
	   string twenty;

	   int user_count = 0;


	   while (data.good() && user_count<user_total)
	   {

		   getline(data, one, ',');
		   getline(data, two, ',');
		   getline(data, three, ',');
		   getline(data, four, ',');
		   getline(data, five, ',');
		   getline(data, six, ',');
		   getline(data, seven, ',');
		   getline(data, eight, ',');
		   getline(data, nine, ',');
		   getline(data, ten, ',');
		   getline(data, eleven, ',');
		   getline(data, twelve, ',');
		   getline(data, thirteen, ',');
		   getline(data, fourteen, ',');
		   getline(data, fifteen, ',');
		   getline(data, sixteen, ',');
		   getline(data, seventeen, ',');
		   getline(data, eighteen, ',');
		   getline(data, nineteen, ',');
		   getline(data, twenty, '\n');


		   /*  ---------------------------   */
		   stringstream geek(one);
		   double md1 = 0;
		   geek >> md1;

		   stringstream geek2(two);
		   double md2 = 0;
		   geek2 >> md2;


		   stringstream geek3(three);
		   double md3 = 0;
		   geek3 >> md3;

		   stringstream geek4(four);
		   double md4 = 0;
		   geek4 >> md4;

		   stringstream geek5(five);
		   double md5 = 0;
		   geek5 >> md5;

		   stringstream geek6(six);
		   double md6 = 0;
		   geek6 >> md6;

		   stringstream geek7(seven);
		   double md7 = 0;
		   geek7 >> md7;

		   stringstream geek8(eight);
		   double md8 = 0;
		   geek8 >> md8;

		   stringstream geek9(nine);
		   double md9 = 0;
		   geek9 >> md9;

		   stringstream geek10(ten);
		   double md10 = 0;
		   geek10 >> md10;

		   stringstream geek11(eleven);
		   double md11 = 0;
		   geek11 >> md11;

		   stringstream geek12(twelve);
		   double md12 = 0;
		   geek12 >> md12;

		   stringstream geek13(thirteen);
		   double md13 = 0;
		   geek13 >> md13;

		   stringstream geek14(fourteen);
		   double md14 = 0;
		   geek14 >> md14;

		   stringstream geek15(fifteen);
		   double md15 = 0;
		   geek15 >> md15;

		   stringstream geek16(sixteen);
		   double md16 = 0;
		   geek16 >> md16;

		   stringstream geek17(seventeen);
		   double md17 = 0;
		   geek17 >> md17;

		   stringstream geek18(eighteen);
		   double md18 = 0;
		   geek18 >> md18;

		   stringstream geek19(nineteen);
		   double md19 = 0;
		   geek19 >> md19;

		   stringstream geek20(twenty);
		   double md20 = 0;
		   geek20 >> md20;

		   if (t == 1)
		   {
			   First_Predict[user_count][t].md = md2*100;
			   Second_Predict[user_count][t].md = md3*100;
			   Third_Predict[user_count][t].md = md4*100;
			   Fourth_Predict[user_count][t].md = md5*100;
			   Fifth_Predict[user_count][t].md = md6*100;
		   }
		   if (t == 2)
		   {
			   First_Predict[user_count][t].md = md3*100;
			   Second_Predict[user_count][t].md = md4*100;
			   Third_Predict[user_count][t].md = md5*100;
			   Fourth_Predict[user_count][t].md = md6*100;
			   Fifth_Predict[user_count][t].md = md7*100;
		   }
		   if (t == 3)
		   {
			   First_Predict[user_count][t].md = md4*100;
			   Second_Predict[user_count][t].md = md5*100;
			   Third_Predict[user_count][t].md = md6*100;
			   Fourth_Predict[user_count][t].md = md7*100;
			   Fifth_Predict[user_count][t].md = md8*100;
		   }
		   if (t == 4)
		   {
			   First_Predict[user_count][t].md = md5*100;
			   Second_Predict[user_count][t].md = md6*100;
			   Third_Predict[user_count][t].md = md7*100;
			   Fourth_Predict[user_count][t].md = md8*100;
			   Fifth_Predict[user_count][t].md = md9*100;
		   }
		   if (t == 5)
		   {
			   First_Predict[user_count][t].md = md6*100;
			   Second_Predict[user_count][t].md = md7*100;
			   Third_Predict[user_count][t].md = md8*100;
			   Fourth_Predict[user_count][t].md = md9*100;
			   Fifth_Predict[user_count][t].md = md10*100;
		   }
		   if (t == 6)
		   {
			   First_Predict[user_count][t].md = md7*100;
			   Second_Predict[user_count][t].md = md8*100;
			   Third_Predict[user_count][t].md = md9*100;
			   Fourth_Predict[user_count][t].md = md10*100;
			   Fifth_Predict[user_count][t].md = md11*100;
		   }
		   if (t == 7)
		   {
			   First_Predict[user_count][t].md = md8*100;
			   Second_Predict[user_count][t].md = md9*100;
			   Third_Predict[user_count][t].md = md10*100;
			   Fourth_Predict[user_count][t].md = md11*100;
			   Fifth_Predict[user_count][t].md = md12*100;
		   }
		   if (t == 8)
		   {
			   First_Predict[user_count][t].md = md9*100;
			   Second_Predict[user_count][t].md = md10*100;
			   Third_Predict[user_count][t].md = md11*100;
			   Fourth_Predict[user_count][t].md = md12*100;
			   Fifth_Predict[user_count][t].md = md13*100;
		   }
		   if (t == 9)
		   {
			   First_Predict[user_count][t].md = md10*100;
			   Second_Predict[user_count][t].md = md11*100;
			   Third_Predict[user_count][t].md = md12*100;
			   Fourth_Predict[user_count][t].md = md13*100;
			   Fifth_Predict[user_count][t].md = md14*100;
		   }
		   if (t == 10)
		   {
			   First_Predict[user_count][t].md = md11*100;
			   Second_Predict[user_count][t].md = md12*100;
			   Third_Predict[user_count][t].md = md13*100;
			   Fourth_Predict[user_count][t].md = md14*100;
			   Fifth_Predict[user_count][t].md = md15*100;
		   }
		   if (t == 11)
		   {
			   First_Predict[user_count][t].md = md12*100;
			   Second_Predict[user_count][t].md = md13*100;
			   Third_Predict[user_count][t].md = md14*100;
			   Fourth_Predict[user_count][t].md = md15*100;
			   Fifth_Predict[user_count][t].md = md16*100;
		   }
		   if (t == 12)
		   {
			   First_Predict[user_count][t].md = md13*100;
			   Second_Predict[user_count][t].md = md14*100;
			   Third_Predict[user_count][t].md = md15*100;
			   Fourth_Predict[user_count][t].md = md16*100;
			   Fifth_Predict[user_count][t].md = md17*100;
		   }
		   if (t == 13)
		   {
			   First_Predict[user_count][t].md = md14*100;
			   Second_Predict[user_count][t].md = md15*100;
			   Third_Predict[user_count][t].md = md16*100;
			   Fourth_Predict[user_count][t].md = md17*100;
			   Fifth_Predict[user_count][t].md = md18*100;
		   }
		   if (t == 14)
		   {
			   First_Predict[user_count][t].md = md15*100;
			   Second_Predict[user_count][t].md = md16*100;
			   Third_Predict[user_count][t].md = md17*100;
			   Fourth_Predict[user_count][t].md = md18*100;
			   Fifth_Predict[user_count][t].md = md19*100;
		   }
		   if (t == 15)
		   {
			   First_Predict[user_count][t].md = md16*100;
			   Second_Predict[user_count][t].md = md17*100;
			   Third_Predict[user_count][t].md = md18*100;
			   Fourth_Predict[user_count][t].md = md19*100;
			   Fifth_Predict[user_count][t].md = md20*100;
		   }
		   if (t == 16)
		   {
			   First_Predict[user_count][t].md = md17*100;
			   Second_Predict[user_count][t].md = md18*100;
			   Third_Predict[user_count][t].md = md19*100;
			   Fourth_Predict[user_count][t].md = md20*100;
		   }
		   if (t == 17)
		   {
			   First_Predict[user_count][t].md = md18*100;
			   Second_Predict[user_count][t].md = md19*100;
			   Third_Predict[user_count][t].md = md20*100;
		   }
		   if (t == 18)
		   {
			   First_Predict[user_count][t].md = md19*100;
			   Second_Predict[user_count][t].md = md20*100;
		   }
		   if (t == 19)
		   {
			   First_Predict[user_count][t].md = md20*100;
		   }
		   if (t == 20)
		   {
			   //First_Predict[user_count][t].b = b2;
			   //Second_Predict[user_count][t].b = b3;
			   //Third_Predict[user_count][t].b = b4;
			   //Fourth_Predict[user_count][t].b = b5;
			   //Fifth_Predict[user_count][t].b = b6;
		   }

		   user_count++;
	   }
	   t++;
   }

   /* p */
   t = 1;
   
   while (t <= 19)
   {
	   stringstream out;
	   out << t;
	   ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\210\\Results\\p\\results_Sample_Data_210_p_" + out.str() + ".csv");

	   if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	   string one;
	   string two;
	   string three;
	   string four;
	   string five;
	   string six;
	   string seven;
	   string eight;
	   string nine;
	   string ten;
	   string eleven;
	   string twelve;
	   string thirteen;
	   string fourteen;
	   string fifteen;
	   string sixteen;
	   string seventeen;
	   string eighteen;
	   string nineteen;
	   string twenty;

	   int user_count = 0;

	   while (data.good() && user_count<user_total)
	   {

		   getline(data, one, ',');
		   getline(data, two, ',');
		   getline(data, three, ',');
		   getline(data, four, ',');
		   getline(data, five, ',');
		   getline(data, six, ',');
		   getline(data, seven, ',');
		   getline(data, eight, ',');
		   getline(data, nine, ',');
		   getline(data, ten, ',');
		   getline(data, eleven, ',');
		   getline(data, twelve, ',');
		   getline(data, thirteen, ',');
		   getline(data, fourteen, ',');
		   getline(data, fifteen, ',');
		   getline(data, sixteen, ',');
		   getline(data, seventeen, ',');
		   getline(data, eighteen, ',');
		   getline(data, nineteen, ',');
		   getline(data, twenty, '\n');


		   /*  ---------------------------   */
		   stringstream geek(one);
		   double p1 = 0;
		   geek >> p1;

		   stringstream geek2(two);
		   double p2 = 0;
		   geek2 >> p2;


		   stringstream geek3(three);
		   double p3 = 0;
		   geek3 >> p3;

		   stringstream geek4(four);
		   double p4 = 0;
		   geek4 >> p4;

		   stringstream geek5(five);
		   double p5 = 0;
		   geek5 >> p5;

		   stringstream geek6(six);
		   double p6 = 0;
		   geek6 >> p6;

		   stringstream geek7(seven);
		   double p7 = 0;
		   geek7 >> p7;

		   stringstream geek8(eight);
		   double p8 = 0;
		   geek8 >> p8;

		   stringstream geek9(nine);
		   double p9 = 0;
		   geek9 >> p9;

		   stringstream geek10(ten);
		   double p10 = 0;
		   geek10 >> p10;

		   stringstream geek11(eleven);
		   double p11 = 0;
		   geek11 >> p11;

		   stringstream geek12(twelve);
		   double p12 = 0;
		   geek12 >> p12;

		   stringstream geek13(thirteen);
		   double p13 = 0;
		   geek13 >> p13;

		   stringstream geek14(fourteen);
		   double p14 = 0;
		   geek14 >> p14;

		   stringstream geek15(fifteen);
		   double p15 = 0;
		   geek15 >> p15;

		   stringstream geek16(sixteen);
		   double p16 = 0;
		   geek16 >> p16;

		   stringstream geek17(seventeen);
		   double p17 = 0;
		   geek17 >> p17;

		   stringstream geek18(eighteen);
		   double p18 = 0;
		   geek18 >> p18;

		   stringstream geek19(nineteen);
		   double p19 = 0;
		   geek19 >> p19;

		   stringstream geek20(twenty);
		   double p20 = 0;
		   geek20 >> p20;

		   if (t == 1)
		   {
			   First_Predict[user_count][t].p = p2;
			   Second_Predict[user_count][t].p = p3;
			   Third_Predict[user_count][t].p = p4;
			   Fourth_Predict[user_count][t].p = p5;
			   Fifth_Predict[user_count][t].p = p6;
		   }
		   if (t == 2)
		   {
			   First_Predict[user_count][t].p = p3;
			   Second_Predict[user_count][t].p = p4;
			   Third_Predict[user_count][t].p = p5;
			   Fourth_Predict[user_count][t].p = p6;
			   Fifth_Predict[user_count][t].p = p7;
		   }
		   if (t == 3)
		   {
			   First_Predict[user_count][t].p = p4;
			   Second_Predict[user_count][t].p = p5;
			   Third_Predict[user_count][t].p = p6;
			   Fourth_Predict[user_count][t].p = p7;
			   Fifth_Predict[user_count][t].p = p8;
		   }
		   if (t == 4)
		   {
			   First_Predict[user_count][t].p = p5;
			   Second_Predict[user_count][t].p = p6;
			   Third_Predict[user_count][t].p = p7;
			   Fourth_Predict[user_count][t].p = p8;
			   Fifth_Predict[user_count][t].p = p9;
		   }
		   if (t == 5)
		   {
			   First_Predict[user_count][t].p = p6;
			   Second_Predict[user_count][t].p = p7;
			   Third_Predict[user_count][t].p = p8;
			   Fourth_Predict[user_count][t].p = p9;
			   Fifth_Predict[user_count][t].p = p10;
		   }
		   if (t == 6)
		   {
			   First_Predict[user_count][t].p = p7;
			   Second_Predict[user_count][t].p = p8;
			   Third_Predict[user_count][t].p = p9;
			   Fourth_Predict[user_count][t].p = p10;
			   Fifth_Predict[user_count][t].p = p11;
		   }
		   if (t == 7)
		   {
			   First_Predict[user_count][t].p = p8;
			   Second_Predict[user_count][t].p = p9;
			   Third_Predict[user_count][t].p = p10;
			   Fourth_Predict[user_count][t].p = p11;
			   Fifth_Predict[user_count][t].p = p12;
		   }
		   if (t == 8)
		   {
			   First_Predict[user_count][t].p = p9;
			   Second_Predict[user_count][t].p = p10;
			   Third_Predict[user_count][t].p = p11;
			   Fourth_Predict[user_count][t].p = p12;
			   Fifth_Predict[user_count][t].p = p13;
		   }
		   if (t == 9)
		   {
			   First_Predict[user_count][t].p = p10;
			   Second_Predict[user_count][t].p = p11;
			   Third_Predict[user_count][t].p = p12;
			   Fourth_Predict[user_count][t].p = p13;
			   Fifth_Predict[user_count][t].p = p14;
		   }
		   if (t == 10)
		   {
			   First_Predict[user_count][t].p = p11;
			   Second_Predict[user_count][t].p = p12;
			   Third_Predict[user_count][t].p = p13;
			   Fourth_Predict[user_count][t].p = p14;
			   Fifth_Predict[user_count][t].p = p15;
		   }
		   if (t == 11)
		   {
			   First_Predict[user_count][t].p = p12;
			   Second_Predict[user_count][t].p = p13;
			   Third_Predict[user_count][t].p = p14;
			   Fourth_Predict[user_count][t].p = p15;
			   Fifth_Predict[user_count][t].p = p16;
		   }
		   if (t == 12)
		   {
			   First_Predict[user_count][t].p = p13;
			   Second_Predict[user_count][t].p = p14;
			   Third_Predict[user_count][t].p = p15;
			   Fourth_Predict[user_count][t].p = p16;
			   Fifth_Predict[user_count][t].p = p17;
		   }
		   if (t == 13)
		   {
			   First_Predict[user_count][t].p = p14;
			   Second_Predict[user_count][t].p = p15;
			   Third_Predict[user_count][t].p = p16;
			   Fourth_Predict[user_count][t].p = p17;
			   Fifth_Predict[user_count][t].p = p18;
		   }
		   if (t == 14)
		   {
			   First_Predict[user_count][t].p = p15;
			   Second_Predict[user_count][t].p = p16;
			   Third_Predict[user_count][t].p = p17;
			   Fourth_Predict[user_count][t].p = p18;
			   Fifth_Predict[user_count][t].p = p19;
		   }
		   if (t == 15)
		   {
			   First_Predict[user_count][t].p = p16;
			   Second_Predict[user_count][t].p = p17;
			   Third_Predict[user_count][t].p = p18;
			   Fourth_Predict[user_count][t].p = p19;
			   Fifth_Predict[user_count][t].p = p20;
		   }
		   if (t == 16)
		   {
			   First_Predict[user_count][t].p = p17;
			   Second_Predict[user_count][t].p = p18;
			   Third_Predict[user_count][t].p = p19;
			   Fourth_Predict[user_count][t].p = p20;
		   }
		   if (t == 17)
		   {
			   First_Predict[user_count][t].p = p18;
			   Second_Predict[user_count][t].p = p19;
			   Third_Predict[user_count][t].p = p20;
		   }
		   if (t == 18)
		   {
			   First_Predict[user_count][t].p = p19;
			   Second_Predict[user_count][t].p = p20;
		   }
		   if (t == 19)
		   {
			   First_Predict[user_count][t].p = p20;
		   }
		   if (t == 20)
		   {
	
		   }

		   user_count++;
	   }
	   t++;
   }

   /* w */
   t = 1;
   
   while (t <= 19)
   {
	   stringstream out;
	   out << t;
	   ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\210\\Results\\w\\results_Sample_Data_210_w_" + out.str() + ".csv");

	   if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	   string one;
	   string two;
	   string three;
	   string four;
	   string five;
	   string six;
	   string seven;
	   string eight;
	   string nine;
	   string ten;
	   string eleven;
	   string twelve;
	   string thirteen;
	   string fourteen;
	   string fifteen;
	   string sixteen;
	   string seventeen;
	   string eighteen;
	   string nineteen;
	   string twenty;

	   int user_count = 0;


	   while (data.good() && user_count<user_total)
	   {

		   getline(data, one, ',');
		   getline(data, two, ',');
		   getline(data, three, ',');
		   getline(data, four, ',');
		   getline(data, five, ',');
		   getline(data, six, ',');
		   getline(data, seven, ',');
		   getline(data, eight, ',');
		   getline(data, nine, ',');
		   getline(data, ten, ',');
		   getline(data, eleven, ',');
		   getline(data, twelve, ',');
		   getline(data, thirteen, ',');
		   getline(data, fourteen, ',');
		   getline(data, fifteen, ',');
		   getline(data, sixteen, ',');
		   getline(data, seventeen, ',');
		   getline(data, eighteen, ',');
		   getline(data, nineteen, ',');
		   getline(data, twenty, '\n');


		   /*  ---------------------------   */
		   stringstream geek(one);
		   double w1 = 0;
		   geek >> w1;

		   stringstream geek2(two);
		   double w2 = 0;
		   geek2 >> w2;

		   stringstream geek3(three);
		   double w3 = 0;
		   geek3 >> w3;

		   stringstream geek4(four);
		   double w4 = 0;
		   geek4 >> w4;

		   stringstream geek5(five);
		   double w5 = 0;
		   geek5 >> w5;

		   stringstream geek6(six);
		   double w6 = 0;
		   geek6 >> w6;

		   stringstream geek7(seven);
		   double w7 = 0;
		   geek7 >> w7;

		   stringstream geek8(eight);
		   double w8 = 0;
		   geek8 >> w8;

		   stringstream geek9(nine);
		   double w9 = 0;
		   geek9 >> w9;

		   stringstream geek10(ten);
		   double w10 = 0;
		   geek10 >> w10;

		   stringstream geek11(eleven);
		   double w11 = 0;
		   geek11 >> w11;

		   stringstream geek12(twelve);
		   double w12 = 0;
		   geek12 >> w12;

		   stringstream geek13(thirteen);
		   double w13 = 0;
		   geek13 >> w13;

		   stringstream geek14(fourteen);
		   double w14 = 0;
		   geek14 >> w14;

		   stringstream geek15(fifteen);
		   double w15 = 0;
		   geek15 >> w15;

		   stringstream geek16(sixteen);
		   double w16 = 0;
		   geek16 >> w16;

		   stringstream geek17(seventeen);
		   double w17 = 0;
		   geek17 >> w17;

		   stringstream geek18(eighteen);
		   double w18 = 0;
		   geek18 >> w18;

		   stringstream geek19(nineteen);
		   double w19 = 0;
		   geek19 >> w19;

		   stringstream geek20(twenty);
		   double w20 = 0;
		   geek20 >> w20;

		   if (t == 1)
		   {
			   First_Predict[user_count][t].w = w2*100;
			   Second_Predict[user_count][t].w = w3*100;
			   Third_Predict[user_count][t].w = w4*100;
			   Fourth_Predict[user_count][t].w = w5*100;
			   Fifth_Predict[user_count][t].w = w6*100;
		   }
		   if (t == 2)
		   {
			   First_Predict[user_count][t].w = w3*100;
			   Second_Predict[user_count][t].w = w4*100;
			   Third_Predict[user_count][t].w = w5*100;
			   Fourth_Predict[user_count][t].w = w6*100;
			   Fifth_Predict[user_count][t].w = w7*100;
		   }
		   if (t == 3)
		   {
			   First_Predict[user_count][t].w = w4*100;
			   Second_Predict[user_count][t].w = w5*100;
			   Third_Predict[user_count][t].w = w6*100;
			   Fourth_Predict[user_count][t].w = w7*100;
			   Fifth_Predict[user_count][t].w = w8*100;
		   }
		   if (t == 4)
		   {
			   First_Predict[user_count][t].w = w5*100;
			   Second_Predict[user_count][t].w = w6*100;
			   Third_Predict[user_count][t].w = w7*100;
			   Fourth_Predict[user_count][t].w = w8*100;
			   Fifth_Predict[user_count][t].w = w9*100;
		   }
		   if (t == 5)
		   {
			   First_Predict[user_count][t].w = w6*100;
			   Second_Predict[user_count][t].w = w7*100;
			   Third_Predict[user_count][t].w = w8*100;
			   Fourth_Predict[user_count][t].w = w9*100;
			   Fifth_Predict[user_count][t].w = w10*100;
		   }
		   if (t == 6)
		   {
			   First_Predict[user_count][t].w = w7*100;
			   Second_Predict[user_count][t].w = w8*100;
			   Third_Predict[user_count][t].w = w9*100;
			   Fourth_Predict[user_count][t].w = w10*100;
			   Fifth_Predict[user_count][t].w = w11*100;
		   }
		   if (t == 7)
		   {
			   First_Predict[user_count][t].w = w8*100;
			   Second_Predict[user_count][t].w = w9*100;
			   Third_Predict[user_count][t].w = w10*100;
			   Fourth_Predict[user_count][t].w = w11*100;
			   Fifth_Predict[user_count][t].w = w12*100;
		   }
		   if (t == 8)
		   {
			   First_Predict[user_count][t].w = w9*100;
			   Second_Predict[user_count][t].w = w10*100;
			   Third_Predict[user_count][t].w = w11*100;
			   Fourth_Predict[user_count][t].w = w12*100;
			   Fifth_Predict[user_count][t].w = w13*100;
		   }
		   if (t == 9)
		   {
			   First_Predict[user_count][t].w = w10*100;
			   Second_Predict[user_count][t].w = w11*100;
			   Third_Predict[user_count][t].w = w12*100;
			   Fourth_Predict[user_count][t].w = w13*100;
			   Fifth_Predict[user_count][t].w = w14*100;
		   }
		   if (t == 10)
		   {
			   First_Predict[user_count][t].w = w11*100;
			   Second_Predict[user_count][t].w = w12*100;
			   Third_Predict[user_count][t].w = w13*100;
			   Fourth_Predict[user_count][t].w = w14*100;
			   Fifth_Predict[user_count][t].w = w15*100;
		   }
		   if (t == 11)
		   {
			   First_Predict[user_count][t].w = w12*100;
			   Second_Predict[user_count][t].w = w13*100;
			   Third_Predict[user_count][t].w = w14*100;
			   Fourth_Predict[user_count][t].w = w15*100;
			   Fifth_Predict[user_count][t].w = w16*100;
		   }
		   if (t == 12)
		   {
			   First_Predict[user_count][t].w = w13*100;
			   Second_Predict[user_count][t].w = w14*100;
			   Third_Predict[user_count][t].w = w15*100;
			   Fourth_Predict[user_count][t].w = w16*100;
			   Fifth_Predict[user_count][t].w = w17*100;
		   }
		   if (t == 13)
		   {
			   First_Predict[user_count][t].w = w14*100;
			   Second_Predict[user_count][t].w = w15*100;
			   Third_Predict[user_count][t].w = w16*100;
			   Fourth_Predict[user_count][t].w = w17*100;
			   Fifth_Predict[user_count][t].w = w18*100;
		   }
		   if (t == 14)
		   {
			   First_Predict[user_count][t].w = w15*100;
			   Second_Predict[user_count][t].w = w16*100;
			   Third_Predict[user_count][t].w = w17*100;
			   Fourth_Predict[user_count][t].w = w18*100;
			   Fifth_Predict[user_count][t].w = w19*100;
		   }
		   if (t == 15)
		   {
			   First_Predict[user_count][t].w = w16*100;
			   Second_Predict[user_count][t].w = w17*100;
			   Third_Predict[user_count][t].w = w18*100;
			   Fourth_Predict[user_count][t].w = w19*100;
			   Fifth_Predict[user_count][t].w = w20*100;
		   }
		   if (t == 16)
		   {
			   First_Predict[user_count][t].w = w17*100;
			   Second_Predict[user_count][t].w = w18*100;
			   Third_Predict[user_count][t].w = w19*100;
			   Fourth_Predict[user_count][t].w = w20*100;
		   }
		   if (t == 17)
		   {
			   First_Predict[user_count][t].w = w18*100;
			   Second_Predict[user_count][t].w = w19*100;
			   Third_Predict[user_count][t].w = w20*100;
		   }
		   if (t == 18)
		   {
			   First_Predict[user_count][t].w = w19*100;
			   Second_Predict[user_count][t].w = w20*100;
		   }
		   if (t == 19)
		   {
			   First_Predict[user_count][t].w = w20*100;

		   }
		   if (t == 20)
		   {
			
		   }

		   user_count++;
	   }
	   t++;
   }

   /* x */
   t = 1;

   while (t <= 19)
   {
	   stringstream out;
	   out << t;
	   ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\210\\Results\\x\\results_Sample_Data_210_x_" + out.str() + ".csv");

	   if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	   string one;
	   string two;
	   string three;
	   string four;
	   string five;
	   string six;
	   string seven;
	   string eight;
	   string nine;
	   string ten;
	   string eleven;
	   string twelve;
	   string thirteen;
	   string fourteen;
	   string fifteen;
	   string sixteen;
	   string seventeen;
	   string eighteen;
	   string nineteen;
	   string twenty;

	   int user_count = 0;


	   while (data.good() && user_count<user_total)
	   {

		   getline(data, one, ',');
		   getline(data, two, ',');
		   getline(data, three, ',');
		   getline(data, four, ',');
		   getline(data, five, ',');
		   getline(data, six, ',');
		   getline(data, seven, ',');
		   getline(data, eight, ',');
		   getline(data, nine, ',');
		   getline(data, ten, ',');
		   getline(data, eleven, ',');
		   getline(data, twelve, ',');
		   getline(data, thirteen, ',');
		   getline(data, fourteen, ',');
		   getline(data, fifteen, ',');
		   getline(data, sixteen, ',');
		   getline(data, seventeen, ',');
		   getline(data, eighteen, ',');
		   getline(data, nineteen, ',');
		   getline(data, twenty, '\n');


		   /*  ---------------------------   */
		   stringstream geek(one);
		   double x1 = 0;
		   geek >> x1;

		   stringstream geek2(two);
		   double x2 = 0;
		   geek2 >> x2;


		   stringstream geek3(three);
		   double x3 = 0;
		   geek3 >> x3;

		   stringstream geek4(four);
		   double x4 = 0;
		   geek4 >> x4;

		   stringstream geek5(five);
		   double x5 = 0;
		   geek5 >> x5;

		   stringstream geek6(six);
		   double x6 = 0;
		   geek6 >> x6;

		   stringstream geek7(seven);
		   double x7 = 0;
		   geek7 >> x7;

		   stringstream geek8(eight);
		   double x8 = 0;
		   geek8 >> x8;

		   stringstream geek9(nine);
		   double x9 = 0;
		   geek9 >> x9;

		   stringstream geek10(ten);
		   double x10 = 0;
		   geek10 >> x10;

		   stringstream geek11(eleven);
		   double x11 = 0;
		   geek11 >> x11;

		   stringstream geek12(twelve);
		   double x12 = 0;
		   geek12 >> x12;

		   stringstream geek13(thirteen);
		   double x13 = 0;
		   geek13 >> x13;

		   stringstream geek14(fourteen);
		   double x14 = 0;
		   geek14 >> x14;

		   stringstream geek15(fifteen);
		   double x15 = 0;
		   geek15 >> x15;

		   stringstream geek16(sixteen);
		   double x16 = 0;
		   geek16 >> x16;

		   stringstream geek17(seventeen);
		   double x17 = 0;
		   geek17 >> x17;

		   stringstream geek18(eighteen);
		   double x18 = 0;
		   geek18 >> x18;

		   stringstream geek19(nineteen);
		   double x19 = 0;
		   geek19 >> x19;

		   stringstream geek20(twenty);
		   double x20 = 0;
		   geek20 >> x20;

		   if (t == 1)
		   {
			   First_Predict[user_count][t].x = x2;
			   Second_Predict[user_count][t].x = x3;
			   Third_Predict[user_count][t].x = x4;
			   Fourth_Predict[user_count][t].x = x5;
			   Fifth_Predict[user_count][t].x = x6;
		   }
		   if (t == 2)
		   {
			   First_Predict[user_count][t].x = x3;
			   Second_Predict[user_count][t].x = x4;
			   Third_Predict[user_count][t].x = x5;
			   Fourth_Predict[user_count][t].x = x6;
			   Fifth_Predict[user_count][t].x = x7;
		   }
		   if (t == 3)
		   {
			   First_Predict[user_count][t].x = x4;
			   Second_Predict[user_count][t].x = x5;
			   Third_Predict[user_count][t].x = x6;
			   Fourth_Predict[user_count][t].x = x7;
			   Fifth_Predict[user_count][t].x = x8;
		   }
		   if (t == 4)
		   {
			   First_Predict[user_count][t].x = x5;
			   Second_Predict[user_count][t].x = x6;
			   Third_Predict[user_count][t].x = x7;
			   Fourth_Predict[user_count][t].x = x8;
			   Fifth_Predict[user_count][t].x = x9;
		   }
		   if (t == 5)
		   {
			   First_Predict[user_count][t].x = x6;
			   Second_Predict[user_count][t].x = x7;
			   Third_Predict[user_count][t].x = x8;
			   Fourth_Predict[user_count][t].x = x9;
			   Fifth_Predict[user_count][t].x = x10;
		   }
		   if (t == 6)
		   {
			   First_Predict[user_count][t].x = x7;
			   Second_Predict[user_count][t].x = x8;
			   Third_Predict[user_count][t].x = x9;
			   Fourth_Predict[user_count][t].x = x10;
			   Fifth_Predict[user_count][t].x = x11;
		   }
		   if (t == 7)
		   {
			   First_Predict[user_count][t].x = x8;
			   Second_Predict[user_count][t].x = x9;
			   Third_Predict[user_count][t].x = x10;
			   Fourth_Predict[user_count][t].x = x11;
			   Fifth_Predict[user_count][t].x = x12;
		   }
		   if (t == 8)
		   {
			   First_Predict[user_count][t].x = x9;
			   Second_Predict[user_count][t].x = x10;
			   Third_Predict[user_count][t].x = x11;
			   Fourth_Predict[user_count][t].x = x12;
			   Fifth_Predict[user_count][t].x = x13;
		   }
		   if (t == 9)
		   {
			   First_Predict[user_count][t].x = x10;
			   Second_Predict[user_count][t].x = x11;
			   Third_Predict[user_count][t].x = x12;
			   Fourth_Predict[user_count][t].x = x13;
			   Fifth_Predict[user_count][t].x = x14;
		   }
		   if (t == 10)
		   {
			   First_Predict[user_count][t].x = x11;
			   Second_Predict[user_count][t].x = x12;
			   Third_Predict[user_count][t].x = x13;
			   Fourth_Predict[user_count][t].x = x14;
			   Fifth_Predict[user_count][t].x = x15;
		   }
		   if (t == 11)
		   {
			   First_Predict[user_count][t].x = x12;
			   Second_Predict[user_count][t].x = x13;
			   Third_Predict[user_count][t].x = x14;
			   Fourth_Predict[user_count][t].x = x15;
			   Fifth_Predict[user_count][t].x = x16;
		   }
		   if (t == 12)
		   {
			   First_Predict[user_count][t].x = x13;
			   Second_Predict[user_count][t].x = x14;
			   Third_Predict[user_count][t].x = x15;
			   Fourth_Predict[user_count][t].x = x16;
			   Fifth_Predict[user_count][t].x = x17;
		   }
		   if (t == 13)
		   {
			   First_Predict[user_count][t].x = x14;
			   Second_Predict[user_count][t].x = x15;
			   Third_Predict[user_count][t].x = x16;
			   Fourth_Predict[user_count][t].x = x17;
			   Fifth_Predict[user_count][t].x = x18;
		   }
		   if (t == 14)
		   {
			   First_Predict[user_count][t].x = x15;
			   Second_Predict[user_count][t].x = x16;
			   Third_Predict[user_count][t].x = x17;
			   Fourth_Predict[user_count][t].x = x18;
			   Fifth_Predict[user_count][t].x = x19;
		   }
		   if (t == 15)
		   {
			   First_Predict[user_count][t].x = x16;
			   Second_Predict[user_count][t].x = x17;
			   Third_Predict[user_count][t].x = x18;
			   Fourth_Predict[user_count][t].x = x19;
			   Fifth_Predict[user_count][t].x = x20;
		   }
		   if (t == 16)
		   {
			   First_Predict[user_count][t].x = x17;
			   Second_Predict[user_count][t].x = x18;
			   Third_Predict[user_count][t].x = x19;
			   Fourth_Predict[user_count][t].x = x20;
		   }
		   if (t == 17)
		   {
			   First_Predict[user_count][t].x = x18;
			   Second_Predict[user_count][t].x = x19;
			   Third_Predict[user_count][t].x = x20;
		   }
		   if (t == 18)
		   {
			   First_Predict[user_count][t].x = x19;
			   Second_Predict[user_count][t].x = x20;
		   }
		   if (t == 19)
		   {
			   First_Predict[user_count][t].x = x20;
		   }
		   if (t == 20)
		   {
			   
		   }

		   user_count++;
	   }
	   t++;
   }

   /* y */
   t = 1;

   while (t <= 19)
   {
	   stringstream out;
	   out << t;
	   ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\210\\Results\\y\\results_Sample_Data_210_y_" + out.str() + ".csv");

	   if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	   string one;
	   string two;
	   string three;
	   string four;
	   string five;
	   string six;
	   string seven;
	   string eight;
	   string nine;
	   string ten;
	   string eleven;
	   string twelve;
	   string thirteen;
	   string fourteen;
	   string fifteen;
	   string sixteen;
	   string seventeen;
	   string eighteen;
	   string nineteen;
	   string twenty;

	   int user_count = 0;


	   while (data.good() && user_count<user_total)
	   {

		   getline(data, one, ',');
		   getline(data, two, ',');
		   getline(data, three, ',');
		   getline(data, four, ',');
		   getline(data, five, ',');
		   getline(data, six, ',');
		   getline(data, seven, ',');
		   getline(data, eight, ',');
		   getline(data, nine, ',');
		   getline(data, ten, ',');
		   getline(data, eleven, ',');
		   getline(data, twelve, ',');
		   getline(data, thirteen, ',');
		   getline(data, fourteen, ',');
		   getline(data, fifteen, ',');
		   getline(data, sixteen, ',');
		   getline(data, seventeen, ',');
		   getline(data, eighteen, ',');
		   getline(data, nineteen, ',');
		   getline(data, twenty, '\n');


		   /*  ---------------------------   */
		   stringstream geek(one);
		   double y1 = 0;
		   geek >> y1;

		   stringstream geek2(two);
		   double y2 = 0;
		   geek2 >> y2;


		   stringstream geek3(three);
		   double y3 = 0;
		   geek3 >> y3;

		   stringstream geek4(four);
		   double y4 = 0;
		   geek4 >> y4;

		   stringstream geek5(five);
		   double y5 = 0;
		   geek5 >> y5;

		   stringstream geek6(six);
		   double y6 = 0;
		   geek6 >> y6;

		   stringstream geek7(seven);
		   double y7 = 0;
		   geek7 >> y7;

		   stringstream geek8(eight);
		   double y8 = 0;
		   geek8 >> y8;

		   stringstream geek9(nine);
		   double y9 = 0;
		   geek9 >> y9;

		   stringstream geek10(ten);
		   double y10 = 0;
		   geek10 >> y10;

		   stringstream geek11(eleven);
		   double y11 = 0;
		   geek11 >> y11;

		   stringstream geek12(twelve);
		   double y12 = 0;
		   geek12 >> y12;

		   stringstream geek13(thirteen);
		   double y13 = 0;
		   geek13 >> y13;

		   stringstream geek14(fourteen);
		   double y14 = 0;
		   geek14 >> y14;

		   stringstream geek15(fifteen);
		   double y15 = 0;
		   geek15 >> y15;

		   stringstream geek16(sixteen);
		   double y16 = 0;
		   geek16 >> y16;

		   stringstream geek17(seventeen);
		   double y17 = 0;
		   geek17 >> y17;

		   stringstream geek18(eighteen);
		   double y18 = 0;
		   geek18 >> y18;

		   stringstream geek19(nineteen);
		   double y19 = 0;
		   geek19 >> y19;

		   stringstream geek20(twenty);
		   double y20 = 0;
		   geek20 >> y20;

		   if (t == 1)
		   {
			   First_Predict[user_count][t].y = y2;
			   Second_Predict[user_count][t].y = y3;
			   Third_Predict[user_count][t].y = y4;
			   Fourth_Predict[user_count][t].y = y5;
			   Fifth_Predict[user_count][t].y = y6;
		   }
		   if (t == 2)
		   {
			   First_Predict[user_count][t].y = y3;
			   Second_Predict[user_count][t].y = y4;
			   Third_Predict[user_count][t].y = y5;
			   Fourth_Predict[user_count][t].y = y6;
			   Fifth_Predict[user_count][t].y = y7;
		   }
		   if (t == 3)
		   {
			   First_Predict[user_count][t].y = y4;
			   Second_Predict[user_count][t].y = y5;
			   Third_Predict[user_count][t].y = y6;
			   Fourth_Predict[user_count][t].y = y7;
			   Fifth_Predict[user_count][t].y = y8;
		   }
		   if (t == 4)
		   {
			   First_Predict[user_count][t].y = y5;
			   Second_Predict[user_count][t].y = y6;
			   Third_Predict[user_count][t].y = y7;
			   Fourth_Predict[user_count][t].y = y8;
			   Fifth_Predict[user_count][t].y = y9;
		   }
		   if (t == 5)
		   {
			   First_Predict[user_count][t].y = y6;
			   Second_Predict[user_count][t].y = y7;
			   Third_Predict[user_count][t].y = y8;
			   Fourth_Predict[user_count][t].y = y9;
			   Fifth_Predict[user_count][t].y = y10;
		   }
		   if (t == 6)
		   {
			   First_Predict[user_count][t].y = y7;
			   Second_Predict[user_count][t].y = y8;
			   Third_Predict[user_count][t].y = y9;
			   Fourth_Predict[user_count][t].y = y10;
			   Fifth_Predict[user_count][t].y = y11;
		   }
		   if (t == 7)
		   {
			   First_Predict[user_count][t].y = y8;
			   Second_Predict[user_count][t].y = y9;
			   Third_Predict[user_count][t].y = y10;
			   Fourth_Predict[user_count][t].y = y11;
			   Fifth_Predict[user_count][t].y = y12;
		   }
		   if (t == 8)
		   {
			   First_Predict[user_count][t].y = y9;
			   Second_Predict[user_count][t].y = y10;
			   Third_Predict[user_count][t].y = y11;
			   Fourth_Predict[user_count][t].y = y12;
			   Fifth_Predict[user_count][t].y = y13;
		   }
		   if (t == 9)
		   {
			   First_Predict[user_count][t].y = y10;
			   Second_Predict[user_count][t].y = y11;
			   Third_Predict[user_count][t].y = y12;
			   Fourth_Predict[user_count][t].y = y13;
			   Fifth_Predict[user_count][t].y = y14;
		   }
		   if (t == 10)
		   {
			   First_Predict[user_count][t].y = y11;
			   Second_Predict[user_count][t].y = y12;
			   Third_Predict[user_count][t].y = y13;
			   Fourth_Predict[user_count][t].y = y14;
			   Fifth_Predict[user_count][t].y = y15;
		   }
		   if (t == 11)
		   {
			   First_Predict[user_count][t].y = y12;
			   Second_Predict[user_count][t].y = y13;
			   Third_Predict[user_count][t].y = y14;
			   Fourth_Predict[user_count][t].y = y15;
			   Fifth_Predict[user_count][t].y = y16;
		   }
		   if (t == 12)
		   {
			   First_Predict[user_count][t].y = y13;
			   Second_Predict[user_count][t].y = y14;
			   Third_Predict[user_count][t].y = y15;
			   Fourth_Predict[user_count][t].y = y16;
			   Fifth_Predict[user_count][t].y = y17;
		   }
		   if (t == 13)
		   {
			   First_Predict[user_count][t].y = y14;
			   Second_Predict[user_count][t].y = y15;
			   Third_Predict[user_count][t].y = y16;
			   Fourth_Predict[user_count][t].y = y17;
			   Fifth_Predict[user_count][t].y = y18;
		   }
		   if (t == 14)
		   {
			   First_Predict[user_count][t].y = y15;
			   Second_Predict[user_count][t].y = y16;
			   Third_Predict[user_count][t].y = y17;
			   Fourth_Predict[user_count][t].y = y18;
			   Fifth_Predict[user_count][t].y = y19;
		   }
		   if (t == 15)
		   {
			   First_Predict[user_count][t].y = y16;
			   Second_Predict[user_count][t].y = y17;
			   Third_Predict[user_count][t].y = y18;
			   Fourth_Predict[user_count][t].y = y19;
			   Fifth_Predict[user_count][t].y = y20;
		   }
		   if (t == 16)
		   {
			   First_Predict[user_count][t].y = y17;
			   Second_Predict[user_count][t].y = y18;
			   Third_Predict[user_count][t].y = y19;
			   Fourth_Predict[user_count][t].y = y20;
		   }
		   if (t == 17)
		   {
			   First_Predict[user_count][t].y = y18;
			   Second_Predict[user_count][t].y = y19;
			   Third_Predict[user_count][t].y = y20;
		   }
		   if (t == 18)
		   {
			   First_Predict[user_count][t].y = y19;
			   Second_Predict[user_count][t].y = y20;
		   }
		   if (t == 19)
		   {
			   First_Predict[user_count][t].y = y20;
		   }
		   if (t == 20)
		   {
			  
		   }

		   user_count++;
	   }
	   t++;
   }
   
}

void ReadReal()
{
	/* new dataset --- beginning-----  */
	for (int i = 0;i<user_total;i++)
	{
		users[i] = new Users[200];
		R[i] = new int[200];
		count_time[i] = 200;

	}
  {
	ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\210\\Sample_Data_210+21_b.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;

	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty, '\n');


		/*  ---------------------------   */
		stringstream geek(one);
		double b1 = 0;
		geek >> b1;
		users[user_count][1].b = b1;

		stringstream geek2(two);
		double b2 = 0;
		geek2 >> b2;
		users[user_count][2].b = b2;

		stringstream geek3(three);
		double b3 = 0;
		geek3 >> b3;
		users[user_count][3].b = b3;

		stringstream geek4(four);
		double b4 = 0;
		geek4 >> b4;
		users[user_count][4].b = b4;

		stringstream geek5(five);
		double b5 = 0;
		geek5 >> b5;
		users[user_count][5].b = b5;

		stringstream geek6(six);
		double b6 = 0;
		geek6 >> b6;
		users[user_count][6].b = b6;

		stringstream geek7(seven);
		double b7 = 0;
		geek7 >> b7;
		users[user_count][7].b = b7;

		stringstream geek8(eight);
		double b8 = 0;
		geek8 >> b8;
		users[user_count][8].b = b8;

		stringstream geek9(nine);
		double b9 = 0;
		geek9 >> b9;
		users[user_count][9].b = b9;

		stringstream geek10(ten);
		double b10 = 0;
		geek10 >> b10;
		users[user_count][10].b = b10;

		stringstream geek11(eleven);
		double b11 = 0;
		geek11 >> b11;
		users[user_count][11].b = b11;

		stringstream geek12(twelve);
		double b12 = 0;
		geek12 >> b12;
		users[user_count][12].b = b12;

		stringstream geek13(thirteen);
		double b13 = 0;
		geek13 >> b13;
		users[user_count][13].b = b13;

		stringstream geek14(fourteen);
		double b14 = 0;
		geek14 >> b14;
		users[user_count][14].b = b14;

		stringstream geek15(fifteen);
		double b15 = 0;
		geek15 >> b15;
		users[user_count][15].b = b15;

		stringstream geek16(sixteen);
		double b16 = 0;
		geek16 >> b16;
		users[user_count][16].b = b16;

		stringstream geek17(seventeen);
		double b17 = 0;
		geek17 >> b17;
		users[user_count][17].b = b17;

		stringstream geek18(eighteen);
		double b18 = 0;
		geek18 >> b18;
		users[user_count][18].b = b18;

		stringstream geek19(nineteen);
		double b19 = 0;
		geek19 >> b19;
		users[user_count][19].b = b19;

		stringstream geek20(twenty);
		double b20 = 0;
		geek20 >> b20;
		users[user_count][20].b = b20;

		user_count++;
	}

   }
	
   {
	ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\210\\Sample_Data_210+21_id.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;

	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty, '\n');


		/*  ---------------------------   */
		stringstream geek(one);
		double id1 = 0;
		geek >> id1;
		users[user_count][1].id = id1*100;

		stringstream geek2(two);
		double id2 = 0;
		geek2 >> id2;
		users[user_count][2].id = id2*100;

		stringstream geek3(three);
		double id3 = 0;
		geek3 >> id3;
		users[user_count][3].id = id3*100;

		stringstream geek4(four);
		double id4 = 0;
		geek4 >> id4;
		users[user_count][4].id = id4*100;

		stringstream geek5(five);
		double id5 = 0;
		geek5 >> id5;
		users[user_count][5].id = id5*100;

		stringstream geek6(six);
		double id6 = 0;
		geek6 >> id6;
		users[user_count][6].id = id6*100;

		stringstream geek7(seven);
		double id7 = 0;
		geek7 >> id7;
		users[user_count][7].id = id7*100;

		stringstream geek8(eight);
		double id8 = 0;
		geek8 >> id8;
		users[user_count][8].id = id8*100;

		stringstream geek9(nine);
		double id9 = 0;
		geek9 >> id9;
		users[user_count][9].id = id9*100;

		stringstream geek10(ten);
		double id10 = 0;
		geek10 >> id10;
		users[user_count][10].id = id10*100;

		stringstream geek11(eleven);
		double id11 = 0;
		geek11 >> id11;
		users[user_count][11].id = id11*100;

		stringstream geek12(twelve);
		double id12 = 0;
		geek12 >> id12;
		users[user_count][12].id = id12*100;

		stringstream geek13(thirteen);
		double id13 = 0;
		geek13 >> id13;
		users[user_count][13].id = id13*100;

		stringstream geek14(fourteen);
		double id14 = 0;
		geek14 >> id14;
		users[user_count][14].id = id14*100;

		stringstream geek15(fifteen);
		double id15 = 0;
		geek15 >> id15;
		users[user_count][15].id = id15*100;

		stringstream geek16(sixteen);
		double id16 = 0;
		geek16 >> id16;
		users[user_count][16].id = id16*100;

		stringstream geek17(seventeen);
		double id17 = 0;
		geek17 >> id17;
		users[user_count][17].id = id17*100;

		stringstream geek18(eighteen);
		double id18 = 0;
		geek18 >> id18;
		users[user_count][18].id = id18*100;

		stringstream geek19(nineteen);
		double id19 = 0;
		geek19 >> id19;
		users[user_count][19].id = id19*100;

		stringstream geek20(twenty);
		double id20 = 0;
		geek20 >> id20;
		users[user_count][20].id = id20*100;

		user_count++;
	}
	}

	{
	ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\210\\Sample_Data_210+21_md.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;

	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty, '\n');


		/*  ---------------------------   */
		stringstream geek(one);
		double md1 = 0;
		geek >> md1;
		users[user_count][1].md = md1*100;

		stringstream geek2(two);
		double md2 = 0;
		geek2 >> md2;
		users[user_count][2].md = md2*100;

		stringstream geek3(three);
		double md3 = 0;
		geek3 >> md3;
		users[user_count][3].md = md3*100;

		stringstream geek4(four);
		double md4 = 0;
		geek4 >> md4;
		users[user_count][4].md = md4*100;

		stringstream geek5(five);
		double md5 = 0;
		geek5 >> md5;
		users[user_count][5].md = md5*100;

		stringstream geek6(six);
		double md6 = 0;
		geek6 >> md6;
		users[user_count][6].md = md6*100;

		stringstream geek7(seven);
		double md7 = 0;
		geek7 >> md7;
		users[user_count][7].md = md7*100;

		stringstream geek8(eight);
		double md8 = 0;
		geek8 >> md8;
		users[user_count][8].md = md8*100;

		stringstream geek9(nine);
		double md9 = 0;
		geek9 >> md9;
		users[user_count][9].md = md9*100;

		stringstream geek10(ten);
		double md10 = 0;
		geek10 >> md10;
		users[user_count][10].md = md10*100;

		stringstream geek11(eleven);
		double md11 = 0;
		geek11 >> md11;
		users[user_count][11].md = md11*100;

		stringstream geek12(twelve);
		double md12 = 0;
		geek12 >> md12;
		users[user_count][12].md = md12*100;

		stringstream geek13(thirteen);
		double md13 = 0;
		geek13 >> md13;
		users[user_count][13].md = md13*100;

		stringstream geek14(fourteen);
		double md14 = 0;
		geek14 >> md14;
		users[user_count][14].md = md14*100;

		stringstream geek15(fifteen);
		double md15 = 0;
		geek15 >> md15;
		users[user_count][15].md = md15*100;

		stringstream geek16(sixteen);
		double md16 = 0;
		geek16 >> md16;
		users[user_count][16].md = md16*100;

		stringstream geek17(seventeen);
		double md17 = 0;
		geek17 >> md17;
		users[user_count][17].md = md17*100;

		stringstream geek18(eighteen);
		double md18 = 0;
		geek18 >> md18;
		users[user_count][18].md = md18*100;

		stringstream geek19(nineteen);
		double md19 = 0;
		geek19 >> md19;
		users[user_count][19].md = md19*100;

		stringstream geek20(twenty);
		double md20 = 0;
		geek20 >> md20;
		users[user_count][20].md = md20*100;

		user_count++;
	}
	}

	{
	ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\210\\Sample_Data_210+21_w.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;

	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty, '\n');


		/*  ---------------------------   */
		stringstream geek(one);
		double w1 = 0;
		geek >> w1;
		users[user_count][1].w = w1*100;

		stringstream geek2(two);
		double w2 = 0;
		geek2 >> w2;
		users[user_count][2].w = w2*100;

		stringstream geek3(three);
		double w3 = 0;
		geek3 >> w3;
		users[user_count][3].w = w3*100;

		stringstream geek4(four);
		double w4 = 0;
		geek4 >> w4;
		users[user_count][4].w = w4*100;

		stringstream geek5(five);
		double w5 = 0;
		geek5 >> w5;
		users[user_count][5].w = w5*100;

		stringstream geek6(six);
		double w6 = 0;
		geek6 >> w6;
		users[user_count][6].w = w6*100;

		stringstream geek7(seven);
		double w7 = 0;
		geek7 >> w7;
		users[user_count][7].w = w7*100;

		stringstream geek8(eight);
		double w8 = 0;
		geek8 >> w8;
		users[user_count][8].w = w8*100;

		stringstream geek9(nine);
		double w9 = 0;
		geek9 >> w9;
		users[user_count][9].w = w9*100;

		stringstream geek10(ten);
		double w10 = 0;
		geek10 >> w10;
		users[user_count][10].w = w10*100;

		stringstream geek11(eleven);
		double w11 = 0;
		geek11 >> w11;
		users[user_count][11].w = w11*100;

		stringstream geek12(twelve);
		double w12 = 0;
		geek12 >> w12;
		users[user_count][12].w = w12*100;

		stringstream geek13(thirteen);
		double w13 = 0;
		geek13 >> w13;
		users[user_count][13].w = w13*100;

		stringstream geek14(fourteen);
		double w14 = 0;
		geek14 >> w14;
		users[user_count][14].w = w14*100;

		stringstream geek15(fifteen);
		double w15 = 0;
		geek15 >> w15;
		users[user_count][15].w = w15*100;

		stringstream geek16(sixteen);
		double w16 = 0;
		geek16 >> w16;
		users[user_count][16].w = w16*100;

		stringstream geek17(seventeen);
		double w17 = 0;
		geek17 >> w17;
		users[user_count][17].w = w17*100;

		stringstream geek18(eighteen);
		double w18 = 0;
		geek18 >> w18;
		users[user_count][18].w = w18*100;

		stringstream geek19(nineteen);
		double w19 = 0;
		geek19 >> w19;
		users[user_count][19].w = w19*100;

		stringstream geek20(twenty);
		double w20 = 0;
		geek20 >> w20;
		users[user_count][20].w = w20*100;

		user_count++;
	}
	}

	{
	ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\210\\Sample_Data_210+21_p.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;

	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty, '\n');


		/*  ---------------------------   */
		stringstream geek(one);
		double p1 = 0;
		geek >> p1;
		users[user_count][1].p = p1;

		stringstream geek2(two);
		double p2 = 0;
		geek2 >> p2;
		users[user_count][2].p = p2;

		stringstream geek3(three);
		double p3 = 0;
		geek3 >> p3;
		users[user_count][3].p = p3;

		stringstream geek4(four);
		double p4 = 0;
		geek4 >> p4;
		users[user_count][4].p = p4;

		stringstream geek5(five);
		double p5 = 0;
		geek5 >> p5;
		users[user_count][5].p = p5;

		stringstream geek6(six);
		double p6 = 0;
		geek6 >> p6;
		users[user_count][6].p = p6;

		stringstream geek7(seven);
		double p7 = 0;
		geek7 >> p7;
		users[user_count][7].p = p7;

		stringstream geek8(eight);
		double p8 = 0;
		geek8 >> p8;
		users[user_count][8].p = p8;

		stringstream geek9(nine);
		double p9 = 0;
		geek9 >> p9;
		users[user_count][9].p = p9;

		stringstream geek10(ten);
		double p10 = 0;
		geek10 >> p10;
		users[user_count][10].p = p10;

		stringstream geek11(eleven);
		double p11 = 0;
		geek11 >> p11;
		users[user_count][11].p = p11;

		stringstream geek12(twelve);
		double p12 = 0;
		geek12 >> p12;
		users[user_count][12].p = p12;

		stringstream geek13(thirteen);
		double p13 = 0;
		geek13 >> p13;
		users[user_count][13].p = p13;

		stringstream geek14(fourteen);
		double p14 = 0;
		geek14 >> p14;
		users[user_count][14].p = p14;

		stringstream geek15(fifteen);
		double p15 = 0;
		geek15 >> p15;
		users[user_count][15].p = p15;

		stringstream geek16(sixteen);
		double p16 = 0;
		geek16 >> p16;
		users[user_count][16].p = p16;

		stringstream geek17(seventeen);
		double p17 = 0;
		geek17 >> p17;
		users[user_count][17].p = p17;

		stringstream geek18(eighteen);
		double p18 = 0;
		geek18 >> p18;
		users[user_count][18].p = p18;

		stringstream geek19(nineteen);
		double p19 = 0;
		geek19 >> p19;
		users[user_count][19].p = p19;

		stringstream geek20(twenty);
		double p20 = 0;
		geek20 >> p20;
		users[user_count][20].p = p20;

		user_count++;
	}
	}

	{
	ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\210\\Sample_Data_210+21_x.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;

	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty, '\n');


		/*  ---------------------------   */
		stringstream geek(one);
		double x1 = 0;
		geek >> x1;
		users[user_count][1].x = x1;
		if (x1 <= min_x)
			min_x = x1;
		if (x1>max_x)
			max_x = x1;

		stringstream geek2(two);
		double x2 = 0;
		geek2 >> x2;
		users[user_count][2].x = x2;
		if (x2 <= min_x)
			min_x = x2;
		if (x2>max_x)
			max_x = x2;
		stringstream geek3(three);
		double x3 = 0;
		geek3 >> x3;
		users[user_count][3].x = x3;
		if (x3 <= min_x)
			min_x = x3;
		if (x3>max_x)
			max_x = x3;

		stringstream geek4(four);
		double x4 = 0;
		geek4 >> x4;
		users[user_count][4].x = x4;
		if (x4 <= min_x)
			min_x = x4;
		if (x4>max_x)
			max_x = x4;

		stringstream geek5(five);
		double x5 = 0;
		geek5 >> x5;
		users[user_count][5].x = x5;
		if (x5 <= min_x)
			min_x = x5;
		if (x5>max_x)
			max_x = x5;
		stringstream geek6(six);
		double x6 = 0;
		geek6 >> x6;
		users[user_count][6].x = x6;
		if (x6 <= min_x)
			min_x = x6;
		if (x6>max_x)
			max_x = x6;
		stringstream geek7(seven);
		double x7 = 0;
		geek7 >> x7;
		users[user_count][7].x = x7;
		if (x7 <= min_x)
			min_x = x7;
		if (x7>max_x)
			max_x = x7;
		stringstream geek8(eight);
		double x8 = 0;
		geek8 >> x8;
		users[user_count][8].x = x8;
		if (x8 <= min_x)
			min_x = x8;
		if (x8>max_x)
			max_x = x8;
		stringstream geek9(nine);
		double x9 = 0;
		geek9 >> x9;
		users[user_count][9].x = x9;
		if (x9 <= min_x)
			min_x = x9;
		if (x9>max_x)
			max_x = x9;
		stringstream geek10(ten);
		double x10 = 0;
		geek10 >> x10;
		users[user_count][10].x = x10;
		if (x10 <= min_x)
			min_x = x10;
		if (x10>max_x)
			max_x = x10;
		stringstream geek11(eleven);
		double x11 = 0;
		geek11 >> x11;
		users[user_count][11].x = x11;
		if (x11 <= min_x)
			min_x = x11;
		if (x11>max_x)
			max_x = x11;
		stringstream geek12(twelve);
		double x12 = 0;
		geek12 >> x12;
		users[user_count][12].x = x12;
		if (x12 <= min_x)
			min_x = x12;
		if (x12>max_x)
			max_x = x12;
		stringstream geek13(thirteen);
		double x13 = 0;
		geek13 >> x13;
		users[user_count][13].x = x13;
		if (x13 <= min_x)
			min_x = x13;
		if (x13>max_x)
			max_x = x13;
		stringstream geek14(fourteen);
		double x14 = 0;
		geek14 >> x14;
		users[user_count][14].x = x14;
		if (x14 <= min_x)
			min_x = x14;
		if (x14>max_x)
			max_x = x14;
		stringstream geek15(fifteen);
		double x15 = 0;
		geek15 >> x15;
		users[user_count][15].x = x15;
		if (x15 <= min_x)
			min_x = x15;
		if (x15>max_x)
			max_x = x15;
		stringstream geek16(sixteen);
		double x16 = 0;
		geek16 >> x16;
		users[user_count][16].x = x16;
		if (x16 <= min_x)
			min_x = x16;
		if (x16>max_x)
			max_x = x16;
		stringstream geek17(seventeen);
		double x17 = 0;
		geek17 >> x17;
		users[user_count][17].x = x17;
		if (x17 <= min_x)
			min_x = x17;
		if (x17>max_x)
			max_x = x17;
		stringstream geek18(eighteen);
		double x18 = 0;
		geek18 >> x18;
		users[user_count][18].x = x18;
		if (x18 <= min_x)
			min_x = x18;
		if (x18>max_x)
			max_x = x18;
		stringstream geek19(nineteen);
		double x19 = 0;
		geek19 >> x19;
		users[user_count][19].x = x19;
		if (x19 <= min_x)
			min_x = x19;
		if (x19>max_x)
			max_x = x19;
		stringstream geek20(twenty);
		double x20 = 0;
		geek20 >> x20;
		users[user_count][20].x = x20;
		if (x20 <= min_x)
			min_x = x20;
		if (x20>max_x)
			max_x = x20;
		user_count++;
	}
	}

	{
	ifstream data("C:\\Users\\Erfan\\Desktop\\Results\\210\\Sample_Data_210+21_y.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;

	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty, '\n');


		/*  ---------------------------   */
		stringstream geek(one);
		double y1 = 0;
		geek >> y1;
		users[user_count][1].y = y1;
		if (y1 <= min_y)
			min_y = y1;
		if (y1>max_y)
			max_y = y1;

		stringstream geek2(two);
		double y2 = 0;
		geek2 >> y2;
		users[user_count][2].y = y2;
		if (y2 <= min_y)
			min_y = y2;
		if (y2>max_y)
			max_y = y2;

		stringstream geek3(three);
		double y3 = 0;
		geek3 >> y3;
		users[user_count][3].y = y3;
		if (y3 <= min_y)
			min_y = y3;
		if (y3>max_y)
			max_y = y3;

		stringstream geek4(four);
		double y4 = 0;
		geek4 >> y4;
		users[user_count][4].y = y4;
		if (y4 <= min_y)
			min_y = y4;
		if (y4>max_y)
			max_y = y4;
		stringstream geek5(five);
		double y5 = 0;
		geek5 >> y5;
		users[user_count][5].y = y5;
		if (y5 <= min_y)
			min_y = y5;
		if (y5>max_y)
			max_y = y5;
		stringstream geek6(six);
		double y6 = 0;
		geek6 >> y6;
		users[user_count][6].y = y6;
		if (y6 <= min_y)
			min_y = y6;
		if (y6>max_y)
			max_y = y6;
		stringstream geek7(seven);
		double y7 = 0;
		geek7 >> y7;
		users[user_count][7].y = y7;
		if (y7 <= min_y)
			min_y = y7;
		if (y7>max_y)
			max_y = y7;
		stringstream geek8(eight);
		double y8 = 0;
		geek8 >> y8;
		users[user_count][8].y = y8;
		if (y8 <= min_y)
			min_y = y8;
		if (y8>max_y)
			max_y = y8;
		stringstream geek9(nine);
		double y9 = 0;
		geek9 >> y9;
		users[user_count][9].y = y9;
		if (y9 <= min_y)
			min_y = y9;
		if (y9>max_y)
			max_y = y9;
		stringstream geek10(ten);
		double y10 = 0;
		geek10 >> y10;
		users[user_count][10].y = y10;
		if (y10 <= min_y)
			min_y = y10;
		if (y10>max_y)
			max_y = y10;
		stringstream geek11(eleven);
		double y11 = 0;
		geek11 >> y11;
		users[user_count][11].y = y11;
		if (y11 <= min_y)
			min_y = y11;
		if (y11>max_y)
			max_y = y11;
		stringstream geek12(twelve);
		double y12 = 0;
		geek12 >> y12;
		users[user_count][12].y = y12;
		if (y12 <= min_y)
			min_y = y12;
		if (y12>max_y)
			max_y = y12;
		stringstream geek13(thirteen);
		double y13 = 0;
		geek13 >> y13;
		users[user_count][13].y = y13;
		if (y13 <= min_y)
			min_y = y13;
		if (y13>max_y)
			max_y = y13;
		stringstream geek14(fourteen);
		double y14 = 0;
		geek14 >> y14;
		users[user_count][14].y = y14;
		if (y14 <= min_y)
			min_y = y14;
		if (y14>max_y)
			max_y = y14;
		stringstream geek15(fifteen);
		double y15 = 0;
		geek15 >> y15;
		users[user_count][15].y = y15;
		if (y15 <= min_y)
			min_y = y15;
		if (y15>max_y)
			max_y = y15;
		stringstream geek16(sixteen);
		double y16 = 0;
		geek16 >> y16;
		users[user_count][16].y = y16;
		if (y16 <= min_y)
			min_y = y16;
		if (y16>max_y)
			max_y = y16;
		stringstream geek17(seventeen);
		double y17 = 0;
		geek17 >> y17;
		users[user_count][17].y = y17;
		if (y17 <= min_y)
			min_y = y17;
		if (y17>max_y)
			max_y = y17;
		stringstream geek18(eighteen);
		double y18 = 0;
		geek18 >> y18;
		users[user_count][18].y = y18;
		if (y18 <= min_y)
			min_y = y18;
		if (y18>max_y)
			max_y = y18;
		stringstream geek19(nineteen);
		double y19 = 0;
		geek19 >> y19;
		users[user_count][19].y = y19;
		if (y19 <= min_y)
			min_y = y19;
		if (y19>max_y)
			max_y = y19;
		stringstream geek20(twenty);
		double y20 = 0;
		geek20 >> y20;
		users[user_count][20].y = y20;
		if (y20 <= min_y)
			min_y = y20;
		if (y20>max_y)
			max_y = y20;
		user_count++;
	}
	}

	   
}

void WriteCloudlets()
{
	/* write to files */
	/*------------------------------------------*/
	ofstream myfile;
	myfile.open("C:\\Users\\Erfan\\Desktop\\Cloudlets.csv");
	myfile.precision(15);
	for (int i = 0;i < cloudlet_total;i++)
	{
			myfile << cloudlets[i].x << ",";
			myfile << cloudlets[i].y << ",";
			myfile << cloudlets[i].p << ",";
			myfile << cloudlets[i].b << ",";
		    myfile << endl;
	}

}

void ReadCloudlets()
{
	
		ifstream data("C:\\Users\\Erfan\\Desktop\\Cloudlets.csv");

		if (!data.is_open()) cout << "ERROR: File Open" << '\n';


		string one;
		string two;
		string three;
		string four;
		
		int user_count = 0;

		while (data.good() && user_count<cloudlet_total)
		{

			getline(data, one, ',');
			getline(data, two, ',');
			getline(data, three, ',');
			getline(data, four, '\n');


			/*  ---------------------------   */
			stringstream geek(one);
			double x = 0;
			geek >> x;

			stringstream geek2(two);
			double y = 0;
			geek2 >> y;


			stringstream geek3(three);
			double p = 0;
			geek3 >> p;

			stringstream geek4(four);
			double b = 0;
			geek4 >> b;

			cloudlets[user_count].x = x;
			cloudlets[user_count].y = y;
			cloudlets[user_count].p = p;
			cloudlets[user_count].b = b;

			user_count++;
		}

}
int main()
{
	ReadReal();
    Results_MatrixCompletion();
	
	ReadCloudlets();

	/* allocating memory for V which will be used in OAMC */
	for (int i = 0;i<user_total;i++)
	{
		V[i] = new int[cloudlet_total];
		r[i] = new double[cloudlet_total];
		Temp_V[i] = new int[cloudlet_total];
		F[i] = new int[cloudlet_total];
		N[i] = new int[cloudlet_total];
	}
	/* allocating memory for S which will be used in Assignment */
	for (int i = 0;i<cloudlet_total;i++)
		{
		  S[i] = new bool[user_total];
		  UsageP[i] = new double[21]; /* UsageP[i][t] shows P which has already consumed from cloudlet i at time slot t*/
		  UsageB[i] = new double[21];
		}

    for(int i=0;i<cloudlet_total;i++)
	  for(int t=1;t<=20;t++)
	    {
		 UsageP[i][t] = 0;
		 UsageB[i][t] = 0;
		}
	

	/*beginning of experiments*/
	Initialize();

	auto start3 = high_resolution_clock::now();
	migration = 0;
	cout << OAMC() << " OAMC(WITH PREDICTION AND EVERY W TIME SLOTS AND WITHOUT SAMPLING) IN MICROSECONDS  " << migration << " w= " << w << " load balancing " << load_balance << endl;
	migration = 0;
	load_balance = 0;

	cout<<"Migrations for different datasets\n"<<"New York "<<MigrationNewYork<<" Orlando "<<MigrationOrlando<<" DACT "<<MigrationDACT<<endl;
	MigrationNewYork = MigrationOrlando = MigrationDACT = 0;

	cout << "metrics\n" << endl;

	cout << "computation time " << ComputationTime << endl;
	cout << "offloading time " << OffloadingTime << endl;
	cout << "migration time " << MigrationTime << endl;
	ComputationTime = OffloadingTime = MigrationTime = 0;

	// Get ending time point
	auto stop3 = high_resolution_clock::now();

	auto duration3 = duration_cast<microseconds>(stop3 - start3);

	cout << "Time taken by function: "
		<< duration3.count() << " microseconds" << endl;


	double AveragePConsumption=0;
	double AverageBConsumption=0;

	for (int t = 1;t <= 20;t++)
	{
		for (int j = 0;j<cloudlet_total;j++)
		{
			 AveragePConsumption += (UsageP[j][t])/cloudlets[j].p;
		     AverageBConsumption += (UsageB[j][t]) / cloudlets[j].b;
		}

	}
	
	cout<<"Average p consumption for OAMC: "<<AveragePConsumption/(20*cloudlet_total)<<endl;
	cout << "Average B consumption for OAMC: " << AverageBConsumption / (20 * cloudlet_total) << endl;


	system("pause");
	

	load_balance = 0;
	for (int i = 0;i<cloudlet_total;i++)
		for (int t = 1;t <= 20;t++)
		{
			UsageP[i][t] = 0;
			UsageB[i][t] = 0;
		}

	ComputationTime = 0;
	OffloadingTime = 0;
	MigrationTime = 0;
	Initialize();

	auto start33 = high_resolution_clock::now();
	migration = 0;
	cout << S-OAMC() << " S-OAMC(WITH PREDICTION AND EVERY W TIME SLOTS AND SAMPLING) IN MICROSECONDS  " << migration << " w= " << w << " load balancing " << load_balance << endl;
	migration = 0;
	load_balance = 0;


	cout << "Migrations for different datasets\n" << "New York " << MigrationNewYork << " Orlando " << MigrationOrlando << " DACT " << MigrationDACT << endl;
	MigrationNewYork = MigrationOrlando = MigrationDACT = 0;

	cout << "metrics\n" << endl;

	cout << "computation time " << ComputationTime << endl;
	cout << "offloading time " << OffloadingTime << endl;
	cout << "migration time " << MigrationTime << endl;
	ComputationTime = OffloadingTime = MigrationTime = 0;

	// Get ending time point
	auto stop33 = high_resolution_clock::now();

	auto duration33 = duration_cast<microseconds>(stop33 - start33);

	cout << "Time taken by function: "
		<< duration33.count() << " microseconds" << endl;

	 AveragePConsumption = 0;
	 AverageBConsumption = 0;

	for (int t = 1;t <= 20;t++)
	{
		for (int j = 0;j<cloudlet_total;j++)
		{
			AveragePConsumption += (UsageP[j][t]) / cloudlets[j].p;
			AverageBConsumption += (UsageB[j][t]) / cloudlets[j].b;
		}

	}

	cout << "Average p consumption for S-OAMC: " << AveragePConsumption / (20 * cloudlet_total) << endl;
	cout << "Average B consumption for S-OAMC: " << AverageBConsumption / (20 * cloudlet_total) << endl;

	system("pause");

	load_balance = 0;
	for (int i = 0;i<cloudlet_total;i++)
		for (int t = 1;t <= 20;t++)
		{
			UsageP[i][t] = 0;
			UsageB[i][t] = 0;
		}

	ComputationTime = 0;
	OffloadingTime = 0;
	MigrationTime = 0;

	Initialize();

	auto start44 = high_resolution_clock::now();
	migration1 = 0;
	cout << S-OAMC-WP() << " S-OAMC-WP(WITHOUT PREDICTION AND EVERY TIME SLOT AND WITH SAMPLING)  " << migration1 << " load balancing " << load_balance << endl;

	cout << "Migrations for different datasets\n" << "New York " << MigrationNewYork << " Orlando " << MigrationOrlando << " DACT " << MigrationDACT << endl;
	MigrationNewYork = MigrationOrlando = MigrationDACT = 0;

	cout << "metrics\n" << endl;
	migration1 = 0;

	cout << "computation time " << ComputationTime << endl;
	cout << "offloading time " << OffloadingTime << endl;
	cout << "migration time " << MigrationTime << endl;
	ComputationTime = OffloadingTime = MigrationTime = 0;
	// Get ending time point
	auto stop44 = high_resolution_clock::now();

	auto duration44 = duration_cast<microseconds>(stop44 - start44);

	cout << "Time taken by function: "
		<< duration44.count() << " microseconds" << endl;

	 AveragePConsumption = 0;
	 AverageBConsumption = 0;

	for (int t = 1;t <= 20;t++)
	{
		for (int j = 0;j<cloudlet_total;j++)
		{
			AveragePConsumption += (UsageP[j][t]) / cloudlets[j].p;
			AverageBConsumption += (UsageB[j][t]) / cloudlets[j].b;
		}

	}

	cout << "Average p consumption for S-OAMC-WP: " << AveragePConsumption / (20 * cloudlet_total) << endl;
	cout << "Average B consumption for S-OAMC-WP: " << AverageBConsumption / (20 * cloudlet_total) << endl;

	system("pause");

	load_balance = 0;
	for (int i = 0;i<cloudlet_total;i++)
		for (int t = 1;t <= 20;t++)
		{
			UsageP[i][t] = 0;
			UsageB[i][t] = 0;
		}

	ComputationTime = 0;
	OffloadingTime = 0;
	MigrationTime = 0;

	Initialize();

	auto start99 = high_resolution_clock::now();

	migration = 0;
	cout << G-OAMC() << " G-OAMC(WITH PREDICTION AND EVERY W TIME SLOTS)  " << migration << " load balancing " << load_balance << endl;

	cout << "Migrations for different datasets\n" << "New York " << MigrationNewYork << " Orlando " << MigrationOrlando << " DACT " << MigrationDACT << endl;
	MigrationNewYork = MigrationOrlando = MigrationDACT = 0;

	cout << "metrics\n" << endl;

	cout << "computation time " << ComputationTime << endl;
	cout << "offloading time " << OffloadingTime << endl;
	cout << "migration time " << MigrationTime << endl;
	ComputationTime = OffloadingTime = MigrationTime = 0;
	// Get ending time point
	migration = 0;
	load_balance = 0;
	auto stop99 = high_resolution_clock::now();

	auto duration99 = duration_cast<microseconds>(stop99 - start99);

	cout << "Time taken by function: "
		<< duration99.count() << " microseconds" << endl;

	AveragePConsumption = 0;
	AverageBConsumption = 0;

	for (int t = 1;t <= 20;t++)
	{
		for (int j = 0;j<cloudlet_total;j++)
		{
			AveragePConsumption += (UsageP[j][t]) / cloudlets[j].p;
			AverageBConsumption += (UsageB[j][t]) / cloudlets[j].b;
		}

	}

	cout << "Average p consumption for G-OAMC: " << AveragePConsumption / (20 * cloudlet_total) << endl;
	cout << "Average B consumption for G-OAMC: " << AverageBConsumption / (20 * cloudlet_total) << endl;

	system("pause");


	load_balance = 0;
	for (int i = 0;i<cloudlet_total;i++)
		for (int t = 1;t <= 20;t++)
		{
			UsageP[i][t] = 0;
			UsageB[i][t] = 0;
		}

	ComputationTime = 0;
	OffloadingTime = 0;
	MigrationTime = 0;

	Initialize();

	auto start77 = high_resolution_clock::now();

	migration = 0;
	cout << BFD() << " BFD Migration  " << migration << " load balancing " << load_balance << endl;

	cout << "Migrations for different datasets\n" << "New York " << MigrationNewYork << " Orlando " << MigrationOrlando << " DACT " << MigrationDACT << endl;
	MigrationNewYork = MigrationOrlando = MigrationDACT = 0;

	cout << "metrics\n" << endl;
	migration = 0;
	cout << "computation time " << ComputationTime << endl;
	cout << "offloading time " << OffloadingTime << endl;
	cout << "migration time " << MigrationTime << endl;
	ComputationTime = OffloadingTime = MigrationTime = 0;
	// Get ending time point
	auto stop77 = high_resolution_clock::now();

	auto duration77 = duration_cast<microseconds>(stop77 - start77);

	cout << "Time taken by function: "
		<< duration77.count() << " microseconds" << endl;

	AveragePConsumption = 0;
	AverageBConsumption = 0;

	for (int t = 1;t <= 20;t++)
	{
		for (int j = 0;j<cloudlet_total;j++)
		{
			AveragePConsumption += (UsageP[j][t]) / cloudlets[j].p;
			AverageBConsumption += (UsageB[j][t]) / cloudlets[j].b;
		}

	}

	cout << "Average p consumption for BFD: " << AveragePConsumption / (20 * cloudlet_total) << endl;
	cout << "Average B consumption for BFD: " << AverageBConsumption / (20 * cloudlet_total) << endl;
	system("pause");
	

	load_balance = 0;
	for (int i = 0;i<cloudlet_total;i++)
		for (int t = 1;t <= 20;t++)
		{
			UsageP[i][t] = 0;
			UsageB[i][t] = 0;
		}

	ComputationTime = 0;
	OffloadingTime = 0;
	MigrationTime = 0;

	Initialize();

	
	auto start55 = high_resolution_clock::now();

	IP();

	cout << "Migrations for different datasets\n" << "New York " << MigrationNewYork << " Orlando " << MigrationOrlando << " DACT " << MigrationDACT << endl;
	MigrationNewYork = MigrationOrlando = MigrationDACT = 0;

	cout << "metrics\n" << endl;

	cout << "computation time " << ComputationTime << endl;
	cout << "offloading time " << OffloadingTime << endl;
	cout << "migration time " << MigrationTime << endl;
	ComputationTime = OffloadingTime = MigrationTime = 0;
	// Get ending time point
	auto stop55 = high_resolution_clock::now();

	auto duration55 = duration_cast<microseconds>(stop55 - start55);

	cout << "Time taken by IP: "
		<< duration55.count() << " microseconds" << endl;

	AveragePConsumption = 0;
	AverageBConsumption = 0;

	for (int t = 1;t <= 20;t++)
	{
		for (int j = 0;j<cloudlet_total;j++)
		{
			AveragePConsumption += (UsageP[j][t]) / cloudlets[j].p;
			AverageBConsumption += (UsageB[j][t]) / cloudlets[j].b;
		}

	}

	cout << "Average p consumption for IP: " << AveragePConsumption / (20 * cloudlet_total) << endl;
	cout << "Average B consumption for IP: " << AverageBConsumption / (20 * cloudlet_total) << endl;
	cout << " load balancing for IP " << load_balance << endl;

	system("pause");

	return 0;
}
