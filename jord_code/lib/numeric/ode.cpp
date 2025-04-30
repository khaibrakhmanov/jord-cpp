#include <iostream>
#include <algorithm>
#include <cmath>
#include "ode.h"

TSwitcher::TSwitcher()
{
	value = false;
};

void TSwitcher::TurnOn()
{
	value = true;
};

Tbool TSwitcher::Active()
{
	return value;
}

void TSwitcher::TurnOff()
{
	value = false;
};

//------------------------------------------------

FILE* TBaseODESolver::Get_file_ID()
{
	return pOdeLogFile;
}

void TBaseODESolver::AdjustSolution()
{
	for (Tint i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		sol[i_sol] = sol2[i_sol] + absoluteErrorVector[i_sol];
	};
}

TBaseODESolver::TBaseODESolver()
{
	min_step = 1.0e-10;
	max_step = 0.1;
	solutionFound = true;
	maxIterationsNumber = 20;
}

void TBaseODESolver::Complete()
{
	CloseLog();
};

void TBaseODESolver::ActivateStepSizeControl()
{
	StepSizeControl.TurnOn();
};

Tbool TBaseODESolver::IfStepSizeControl()
{
	return StepSizeControl.Active();
};

Tfloat TBaseODESolver::GetMinStep()
{
	return min_step;
};

Tfloat TBaseODESolver::GetMaxStep()
{
	return max_step;
};

void TBaseODESolver::SetMinStep(const Tfloat arg)
{
	min_step = arg;
};

void TBaseODESolver::SetMaxStep(const Tfloat arg)
{
	max_step = arg;
}

Tbool TBaseODESolver::GetSolutionStatus()
{
	return solutionFound;
}

void TBaseODESolver::CreateLog(char * fname)
{
	//fopen_s(&pOdeLogFile, fname, "w");
	pOdeLogFile = fopen(fname, "w");

	if (!pOdeLogFile)
	{
		std::cout << "unable to open log file for writing!" << std::endl;
	};
}

void TBaseODESolver::UpdateLogSection(char * string)
{
	fprintf(Get_file_ID(), "%s\n\n", string);
}

void TBaseODESolver::CloseLog()
{
	fclose(Get_file_ID());
}

void TBaseODESolver::LogConvergenceParameters()
{
	fprintf(Get_file_ID(), "Convergence parameters at t = %5.2e):\n", GetTime());
	fprintf(Get_file_ID(), "  dt     = %5.2e):\n", GetTimeStep());
	fprintf(Get_file_ID(), "  eps    = %5.2e):\n", GetErrorMax());
	fprintf(Get_file_ID(), "  N_iter = %5.2e):\n", GetIterationsNumber());

	fprintf(Get_file_ID(), "-------------------------------------\n");
};

void TBaseODESolver::SetTime(const Tfloat _time)
{
	time = _time;
};

Tfloat TBaseODESolver::GetTimeStep()
{
	return dt;
};

void TBaseODESolver::DoubleStep()
{
	SetTimeStep(2.0 * GetTimeStep());
};

void TBaseODESolver::HalveStep()
{
	SetTimeStep(0.5 * GetTimeStep());
}

Tfloat TBaseODESolver::GetTime()
{
	return time;
};

void TBaseODESolver::SetTimeStep(const Tfloat _dt)
{
	dt = _dt;
};

void TBaseODESolver::SetTolerance(const Tfloat _tol)
{
	tol = _tol;
};

void TBaseODESolver::ComputeError()
{
	/* iterator pointing to the maximal value in the vector */
	Tfloat_vector::iterator vector_max_Iter;
#ifdef PARALLEL_MODE
#pragma omp parallel num_threads(NUM_THREADS) shared(errorVector, sol, sol2)
	{
#pragma omp for nowait
		for (Tshort i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
		{
			errorVector[i_sol] = fabs((sol[i_sol] - sol2[i_sol]) / (fabs(sol2[i_sol]) + fabs(sol[i_sol])));
		};
	};
#else
	for (Tint i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		//absoluteErrorVector[i_sol] = (sol2[i_sol] - sol[i_sol]) / (pow(2.0, GetOrder()) - 1.0);
		errorVector[i_sol] = fabs((sol2[i_sol] - sol[i_sol]) / std::max(sol[i_sol], sol2[i_sol]));
		absoluteErrorVector[i_sol] = (sol2[i_sol] - sol[i_sol]);
	};
#endif

	vector_max_Iter = std::max_element(errorVector.begin(), errorVector.end());
	SetErrorMax(*vector_max_Iter);
}

void TBaseODESolver::SetErrorMax(const Tfloat err_max)
{
	maximalError = err_max;
}

Tfloat TBaseODESolver::GetErrorMax()
{
	return maximalError;
};

void TBaseODESolver::Configure()
{
	ConfigureDefault();
};

void TBaseODESolver::AllocateSolution()
{
	AllocateVector1D(&initialVector, GetVariablesNumber());
	AllocateVector1D(&sol, GetVariablesNumber());
	AllocateVector1D(&sol2, GetVariablesNumber());
	AllocateVector1D(&solHalf, GetVariablesNumber());
	AllocateVector1D(&errorVector, GetVariablesNumber());
	AllocateVector1D(&absoluteErrorVector, GetVariablesNumber());
};

void TBaseODESolver::CheckSize(const Tint _size)
{
	if (_size != GetVariablesNumber())
		std::cout << "[warning] size of the input vector does not match size of the solution vector" << std::endl;
};

void TBaseODESolver::ReadInitialVector(const Tfloat_vector _inV, const Tint vsize)
{
	CheckSize(vsize);
	initialVector = _inV;
};

Tfloat TBaseODESolver::GetSolution(const Tint index)
{
	return sol[index];
};

Tfloat TBaseODESolver::GetError(const Tint index)
{
	return errorVector[index];
};

Tfloat TBaseODESolver::GetTolerance()
{
	return tol;
};

Tbool TBaseODESolver::Solve()
{
	if (!StepSizeControl.Active())
	{
		/* use fixed step size */
		Step();
	}
	else
	{
		/* use automatic step size control */
		Step();
		Step2();
		ComputeError();
		UserCheck();

		reachedIteration = 0;

		if (GetErrorMax() <= GetTolerance())
		{
			if (GetTimeStep() < GetMaxStep())
			{
				IncreaseNextStep.TurnOn();
				if (richardsonExtrapolation.Active())
				{
					AdjustSolution();
				}
			}
		}
		else
		{
			do
			{
				if (GetTimeStep() < GetMinStep())
				{
					solutionFound = false;
					break;
				}
				else
				{
					SetTimeStep(GetTimeStep() * 0.5);
					Step();
					Step2();
					ComputeError();
					reachedIteration++;
				}
			} while ((GetErrorMax() > GetTolerance()) && (GetIterationsNumber() < GetMaxIterationsNumber()));
			if ((GetIterationsNumber() < GetMaxIterationsNumber()) && richardsonExtrapolation.Active())
				AdjustSolution();
		}

		if (GetIterationsNumber() == GetMaxIterationsNumber())
		{
			LogConvergenceParameters();
		};

	};

	return GetSolutionStatus();
};


Tint TBaseODESolver::GetIterationsNumber()
{
	return reachedIteration;
}

Tint TBaseODESolver::GetMaxIterationsNumber()
{
	return maxIterationsNumber;
}

Tint TBaseODESolver::GetVariablesNumber()
{
	return variablesNumber;
};

void TBaseODESolver::SetVariablesNumber(const Tint arg)
{
	variablesNumber = arg;
};

void TBaseODESolver::ConfigureDefault()
{
	SetVariablesNumber(variablesNumber);
	IncreaseNextStep.TurnOff();
};

void TBaseODESolver::UpdateStep(const Tfloat old_step)
{
	if (IncreaseNextStep.Active())
	{
		SetTimeStep(old_step * 2.0);
		IncreaseNextStep.TurnOff();
	};
};

void TBaseODESolver::ActivateRichardsonExtrapolation()
{
	richardsonExtrapolation.TurnOn();
};

//------------------------------------

TBasicRungeKuttaSolver::TBasicRungeKuttaSolver()
{
	richardsonExtrapolation.TurnOff();
}

Tfloat TBasicRungeKuttaSolver::GetC(const Tint i, const Tint j)
{
	return coeff[i][j];
};

void TBasicRungeKuttaSolver::SetC(const Tint i, const Tint j, const Tfloat arg)
{
	coeff[i][j] = arg;
};

Tfloat_vector_2d* TBasicRungeKuttaSolver::GetCoefficientsPtr()
{
	return &coeff;
}

void TBasicRungeKuttaSolver::AllocateCoefficients()
{
	AllocateVector2D(&coeff, GetOrder(), GetVariablesNumber());
};

Tint TBasicRungeKuttaSolver::GetOrder()
{
	return order;
};

void TBasicRungeKuttaSolver::SetOrder(const Tint arg)
{
	order = arg;
};

//------------- Ordinary RK solver----------------------
TBaseRungeKuttaSolver::TBaseRungeKuttaSolver()
{
	SetOrder(4);
}

void TBaseRungeKuttaSolver::Step()
{
	CalcCoefficients(initialVector, GetTimeStep(), GetTime());

#ifdef PARALLEL_MODE
#pragma omp parallel num_threads(NUM_THREADS) shared(sol, initialVector)
	{
#pragma omp for nowait
		for (int i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
		{
			sol[i_sol] = initialVector[i_sol] + GetTimeStep() / 6.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
		};
	};
#else
	for (Tint i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		sol[i_sol] = initialVector[i_sol] + GetTimeStep() / 6.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
	};
#endif

};

void TBaseRungeKuttaSolver::Step2()
{

	CalcCoefficients(initialVector, 0.5 * GetTimeStep(), GetTime());

#ifdef PARALLEL_MODE
#pragma omp parallel num_threads(NUM_THREADS) shared(solHalf, initialVector)
	{
#pragma omp for nowait
		for (Tshort i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
		{
			solHalf[i_sol] = initialVector[i_sol] + GetTimeStep() / 12.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
		};
	};
#else
	for (Tint i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		solHalf[i_sol] = initialVector[i_sol] + GetTimeStep() / 12.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
	};
#endif

	CalcCoefficients(solHalf, 0.5 * GetTimeStep(), GetTime() + 0.5 * GetTimeStep());

#ifdef PARALLEL_MODE
#pragma omp parallel num_threads(NUM_THREADS) shared(sol2, solHalf)
	{
#pragma omp for nowait
		for (Tshort i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
		{
			sol2[i_sol] = solHalf[i_sol] + GetTimeStep() / 12.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
		};
	};
#else
	for (Tint i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		sol2[i_sol] = solHalf[i_sol] + GetTimeStep() / 12.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
	};
#endif

};
