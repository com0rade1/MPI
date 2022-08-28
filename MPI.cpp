#include <iostream>
#include <mpi.h>
#include <stdlib.h>
#include <thread>
#include <cstdlib>
#include <time.h>

void PairOperationsMethod(int argc, char* argv[])
{
	int proc_num, proc_rank;
	MPI_Status st;
	MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
		MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
		if (proc_rank == 0) {
			int N;
			srand(time(NULL));
			std::cout << "Print N: ";
			double t1 = MPI_Wtime();
			std::cin >> N; 
			int tasksCount = N / proc_num; 
			for (int i = 1; i < proc_num; i++) {
				MPI_Send(&tasksCount, 1, MPI_INT, i, 100, MPI_COMM_WORLD);
			} 
			int* mas = new int[N]; 
			for (int i = 0; i < N; i++) {
				mas[i] = rand() % 1000;
				
			} 

			int temp = 0;
			for (int i = 1; i < proc_num; i++) {
				for (int j = i * tasksCount; j < (i + 1) * tasksCount; j++) {
					temp = mas[j];
					MPI_Send(&temp, 1, MPI_INT, i, 100, MPI_COMM_WORLD);
				}
			}
			int max = mas[0];
			for (int i = 1; i < tasksCount; i++) {
				if (max <= mas[i]) max = mas[i];
			}
			for (int i = 1; i < proc_num; i++) {
				int tempMax = 0;
				MPI_Recv(&tempMax, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD, &st);
				if (max < tempMax) max = tempMax;
			}
			
			double t2 = MPI_Wtime();
			std::cout << "Time: " << t2 - t1;
		}
		else {
			int tasksCount;
			MPI_Recv(&tasksCount, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
			
			int* mas = new int[tasksCount];
			for (int i = 0; i < tasksCount; i++) {
				MPI_Recv(&mas[i], 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
			}
			int max = mas[0];
			for (int i = 0; i < tasksCount; i++) {
				if (max <= mas[i]) max = mas[i];
			}
			MPI_Send(&max, 1, MPI_INT, 0, 100, MPI_COMM_WORLD);
		}
	MPI_Finalize();
	return;
}
void CollectiveOperationsMethod(int argc, char* argv[]) {
	int proc_num, proc_rank, N, tasksCount;
	int* mas = new int[2];int* masProc;
	MPI_Status st;
	MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
		MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
		if (proc_rank == 0) {
			std::cout << "Print N: ";
			std::cin >> N;
			tasksCount = N / proc_num;
		}
		double t1 = MPI_Wtime();
		MPI_Bcast(&tasksCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (proc_rank == 0) {
			mas = new int[N];
			srand(time(NULL));
			for (int i = 0; i < N; i++) {
				mas[i] = rand() % 1000;
				
			}

		}
		masProc = new int[tasksCount];
		MPI_Scatter(&mas[0], tasksCount, MPI_INT, &masProc[0], tasksCount, MPI_INT, 0, MPI_COMM_WORLD);
		int tempMax;
		int max;
		max = masProc[0];
		for (int i = 0; i < tasksCount; i++) {
			if (max < masProc[i]) max = masProc[i];

		}

		MPI_Reduce(&max,&tempMax,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
		double t2 = MPI_Wtime();
		if (proc_rank == 0) {

			std::cout << "Time: " << t2 - t1;
		}
	MPI_Finalize();
	return;
}
void standartMethod(int argc, char* argv[]) {
	std::cout << "Print N: ";
	int N;
	time_t start, end;
	std::cin >> N;
	time(&start);
	int* mas = new int[N];
	srand(time(NULL));
	for (int i = 0; i < N; i++) {
		mas[i] = rand() % 1000;
	}
	int max = mas[0];
	for (int i = 1; i < N; i++) {
		if (max < mas[i]) max = mas[i];
	}
	time(&end);
	double t = difftime(end, start);
	std::cout << "Time: " << t;
}

int main(int argc, char* argv[])
{

	CollectiveOperationsMethod(argc, argv);
}

void CollectiveOperationsMethod(int argc, char* argv[]) {
	int proc_num, proc_rank, N, tasksCount;
	int* mas = new int[2];int* masProc;
	MPI_Status st;
	MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
		MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
		if (proc_rank == 0) {
			std::cout << "Print N: ";
			std::cin >> N;
			tasksCount = N / proc_num;
		}
		double t1 = MPI_Wtime();
		MPI_Bcast(&tasksCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (proc_rank == 0) {
			mas = new int[N];
			srand(time(NULL));
			for (int i = 0; i < N; i++) {
				mas[i] = rand() % 1000;

			}

		}
		masProc = new int[tasksCount];
		MPI_Scatter(&mas[0], tasksCount, MPI_INT, &masProc[0], tasksCount, MPI_INT, 0, MPI_COMM_WORLD);
		int tempMax;
		int max;
		max = masProc[0];
		for (int i = 0; i < tasksCount; i++) {
			if (max < masProc[i]) max = masProc[i];

		}

		MPI_Reduce(&max,&tempMax,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
		double t2 = MPI_Wtime();
		if (proc_rank == 0) {

			std::cout << "Time: " << t2 - t1;
		}
	MPI_Finalize();
	return;
}
void standartMethod(int argc, char* argv[]) {
	std::cout << "Print N: ";
	int N;
	time_t start, end;
	std::cin >> N;
	time(&start);
	int* mas = new int[N];
	srand(time(NULL));
	for (int i = 0; i < N; i++) {
		mas[i] = rand() % 1000;
	}
	int max = mas[0];
	for (int i = 1; i < N; i++) {
		if (max < mas[i]) max = mas[i];
	}
	time(&end);
	double t = difftime(end, start);
	std::cout << "Time: " << t;
}

int main(int argc, char* argv[])
{
	//PairOperationsMethod(argc, argv);
	//standartMethod(argc, argv);
	CollectiveOperationsMethod(argc, argv);
}
