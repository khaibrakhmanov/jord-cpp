#ifndef ODE_INCLUDED
#define ODE_INCLUDED

#include "../coding/utypes.h"

class TSwitcher
{
private:
	Tbool value;
public:
	TSwitcher();

	void TurnOn();
	void TurnOff();
	Tbool Active();
};

class TBaseODESolver
{
private:

	TSwitcher IncreaseNextStep;
	TSwitcher StepSizeControl;

	Tfloat min_step;
	Tfloat max_step;

	Tbool solutionFound;

protected:
	// order of aproximation
	Tshort order;
	TSwitcher richardsonExtrapolation;

	// solution vestors
	Tfloat_vector sol; //obtained using one step
	Tfloat_vector sol2; //obtained using two steps, h/2 and h/2
	Tfloat_vector solHalf; //at the step h / 2
	// initial condition
	Tfloat_vector initialVector;
	// errors vector
	Tfloat_vector errorVector;
	Tfloat_vector absoluteErrorVector;

	Tint variablesNumber;
	// current number of iterations
	Tint reachedIteration;
	Tint maxIterationsNumber;
	// maximal error among the errors for each component of the solution vector
	Tfloat maximalError;

	Tfloat dt;
	Tfloat time;
	// required tolerance
	Tfloat tol;


	//------------- log -----------------
	// Log-file ID
	FILE* pOdeLogFile;
	/* returns ID of the log file */
	FILE* Get_file_ID();
	//-----------------------------------


	void AdjustSolution();

public:
	TBaseODESolver();

	void Complete();

	virtual Tint GetOrder() = 0;
	virtual void SetOrder(const Tint) = 0;

	Tbool Solve();

	virtual void AllocateCoefficients() = 0;
	virtual void Step() = 0; //with h
	virtual void Step2() = 0;// with h/2 + h/2

	void ActivateRichardsonExtrapolation();

	void AllocateSolution();

	void Configure();
	void ConfigureDefault();

	void ActivateStepSizeControl();
	Tbool IfStepSizeControl();

	Tfloat GetTimeStep();
	Tfloat GetTime();
	void SetTime(const Tfloat _time);
	void SetTimeStep(const Tfloat _dt);

	void SetTolerance(const Tfloat);
	Tfloat GetTolerance();
	Tint GetIterationsNumber();
	Tint GetMaxIterationsNumber();

	/* increase step by factor 2 */
	void DoubleStep();
	/* reduce step by factor 2 */
	void HalveStep();
	void UpdateStep(const Tfloat);

	void ReadInitialVector(const Tfloat_vector, const Tint vsize);
	//void ReadInitialVector(const Tfloat);
	void CheckSize(const Tint _size);
	Tfloat GetSolution(const Tint);


	/* returns error for solution with given index */
	Tfloat GetError(const Tint);
	void ComputeError();
	void SetErrorMax(const Tfloat err_max);
	Tfloat GetErrorMax();

	virtual void UserCheck() = 0;

	Tint GetVariablesNumber();
	void SetVariablesNumber(const Tint);

	Tfloat GetMinStep();
	void SetMinStep(const Tfloat);

	Tfloat GetMaxStep();
	void SetMaxStep(const Tfloat);

	/* returns true if the solution is found, false otherwise */
	Tbool GetSolutionStatus();


	//--------- log ------------------

	/* create text file for the errors, messages and warnings reports */
	void CreateLog(char *fname);
	/* write title of the current section to the log file */
	void UpdateLogSection(char *string);
	/* close log file */
	void CloseLog();
	void LogConvergenceParameters();
	//--------------------------------

};

class TBasicRungeKuttaSolver : public TBaseODESolver
{
protected:
	//Tint order;
	Tfloat_vector_2d coeff;

public:
	TBasicRungeKuttaSolver();

	virtual void CalcCoefficients(const Tfloat_vector v0, const Tfloat _dt, const Tfloat _time) = 0;

	virtual void Step() = 0; //with h
	virtual void Step2() = 0;// with h/2 + h/2

	Tfloat GetC(const Tint i, const Tint j);
	void SetC(const Tint i, const Tint j, const Tfloat arg);

	/// <summary>
	/// returns pointer to the array of a scheme coefficients
	/// </summary>
	/// <returns>pointer to coeff</returns>
	Tfloat_vector_2d* GetCoefficientsPtr();

	void AllocateCoefficients();

	Tint GetOrder();
	void SetOrder(const Tint);
};

class TBaseRungeKuttaSolver : public TBasicRungeKuttaSolver
{
public:
	TBaseRungeKuttaSolver();

	virtual void CalcCoefficients(const Tfloat_vector v0, const Tfloat _dt, const Tfloat _time) = 0;

	void Step(); //with h
	void Step2();// with h/2 + h/2
};

#endif
