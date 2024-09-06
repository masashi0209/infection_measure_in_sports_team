#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MT.h"

#define MEMBER 50 // population size
#define FREQ 0.02 // inverse of MEMBER
#define STATES 8 // 0: S, 1: E, 2: P1, 3: P2, 4: Is, 5: Ia, 6: R, 7: quarantined
#define ONE_T 100
#define DELTA 0.01

#define R_0 5.0 // basic reproductive ratio

#define REG_READ_TIME 3

void setRandomSeed()
{
  init_genrand64(1);
}

// rename genrand64_real3 to urand
double urand()
{
  return (double)(genrand64_real3());
}

typedef struct indiv {
	int state; //epidemic states
	int quarantine; //0: in the population, 1: quarantined
	int quarantineDays; // days for quarantine
	int testResult;
	int waitingResult; // 0; not waiting, 1 waiting for result
	int waitingDays; // 0; not waiting, 1 waiting for result
} INDIV;


void infections_in_a_day(int stateNumber[], INDIV indiv[], double beta, double gamma, double rho, double sigma, double eta)
{
	int t, i, member, partner, indivState;
	double rnd, rnd2, force_infection;
	
	
	for(t=0; t<ONE_T; t++){
		force_infection = beta*(double)(stateNumber[2] + stateNumber[3] + stateNumber[4] + stateNumber[5]);
		for(member=0; member<MEMBER; member++){
			indivState = indiv[member].state;
			rnd = urand(); // random real (0, 1)
			switch (indivState) {
			  case 0: //susceptible
				if (rnd < force_infection) { 
					indiv[member].state = 1; 
					stateNumber[0]--;
					stateNumber[1]++;
				}
				break;
			  case 1: // exposed 
				if (rnd < sigma) {//if rnd is smaller than sigma, this individual becomes P1 (state is 2)
					indiv[member].state = 2; 
					stateNumber[1]--;
					stateNumber[2]++;
				}
				
				 break;
			  case 2: // P1
				if (rnd < rho) { //if rnd is smaller than rho, this individual becomes P2
					indiv[member].state = 3; 
					// note that this individuals may already have been quarantined (state 7)
					// change the stateNumber only when this individuals is not quarantined
					if(indiv[member].quarantine == 0){ 
						stateNumber[2]--;
						stateNumber[3]++;
					}
				}
				break;
			  case 3: // P2
				if (rnd < rho) { //P2 individuals will be either Is or Ia
					rnd2 = urand(); //もう一つ乱数を引いて
					if (rnd2 < eta) { 
						indiv[member].state = 4;
						// change the stateNumber only when this individuals is not quarantined
						if(indiv[member].quarantine == 0){
							stateNumber[3]--;
							stateNumber[4]++;
						}
					} else { 
						indiv[member].state = 5;
						// change the stateNumber only when this individuals is not quarantined
						if(indiv[member].quarantine == 0){
							stateNumber[3]--;
							stateNumber[5]++;
						}
					}
				}
				break;
			  case 4: // Is
				if (rnd < gamma) {// recovery?
					indiv[member].state = 6; 
					// change the stateNumber only when this individuals is not quarantined
					if(indiv[member].quarantine == 0){
						stateNumber[4]--;
						stateNumber[6]++;
					}
				}
				break;
			  case 5: // Ia
				if (rnd < gamma) {
					indiv[member].state = 6; 
					// change the stateNumber only when this individuals is not quarantined
					if(indiv[member].quarantine == 0){
						stateNumber[5]--;
						stateNumber[6]++;
					}
				}
				break;
			  case 6: // recovered
				//do nothing
				break;
			} //end switch
		}
		//printf("\n");

	} //one_t
}

void setPCRSensitivity(double PCRSTV[])
{
	PCRSTV[0] = 0;
	PCRSTV[1] = 0;// E
	PCRSTV[2] = 0.33; //P1
	PCRSTV[3] = 0.62; //P2
	PCRSTV[4] = 0.8; //Iss
	PCRSTV[5] = 0.8; //Isa
	PCRSTV[6] = 0; // R
	PCRSTV[7] = 0; // quarantined
}

void setAntigenSensitivity(double antigenSTV[], double PCRSTV[], int scenario)
{
	int i;
	
	if(scenario == 3) for(i=0; i<STATES; i++) antigenSTV[i] = 0.35 * PCRSTV[i];
	else if(scenario == 4) for(i=0; i<STATES; i++) antigenSTV[i] = 0.5 * PCRSTV[i];
	else for(i=0; i<STATES; i++) antigenSTV[i] = 0.70 * PCRSTV[i];
}

void doTest(INDIV indiv[], double sensitivity[])
{
	int member, state;
	double rnd;
	for(member=0; member< MEMBER; member++){
		if(indiv[member].quarantine == 0){ // if the person has not isolated yet, check
			state = indiv[member].state;
			if(sensitivity[state] > 0.0){
				rnd = urand();
				if(rnd < sensitivity[state]){ 
					indiv[member].testResult = 1; // test positive
				}
			}
			// increment testing information
			indiv[member].waitingResult = 1; // they are waiting for test result
			indiv[member].waitingDays = 0; // clear waiting days for test result
		}
	} // member
}

void doAntigenTest(int stateNumber[], INDIV indiv[], double sensitivity[])
{
	int member, state;
	double rnd;
	for(member=0; member< MEMBER; member++){
		if(indiv[member].quarantine == 0){ // if the person has not isolated yet, check
			state = indiv[member].state;
			
			if(sensitivity[state]>0.0){
				rnd = urand();
				if(rnd < sensitivity[state]){
					// reading time is 0, isolate the person immediately
					indiv[member].testResult = 2; // antigen test positive
					indiv[member].quarantine = 1; // isolate
					stateNumber[state]--;
					stateNumber[7]++;
				}
			}
		}
	} // member
}

void initializePopulation(int stateNumber[], INDIV indiv[]){
	int i;

	for(i=0; i<MEMBER; i++){
		indiv[i].state = 0; //all susceptible
		indiv[i].quarantine = 0; //not quarantined
		indiv[i].quarantineDays =0;
		indiv[i].testResult = 0; //0: negative, 1: PCR positive, 2: antigen test positive
		indiv[i].waitingResult = 0; // 0: not waiting, 1: waiting
		indiv[i].waitingDays = 0;
	}
	stateNumber[0] = MEMBER;
	for(i=1; i<STATES; i++) stateNumber[i]=0;
}

void dailySymptomCheck(int stateNumber[], INDIV indiv[]){
	int member;
	
	for(member=0; member< MEMBER; member++){
		if(indiv[member].state == 4){ // this person has a symptom
			if(indiv[member].quarantine == 0){ // if the person has not isolated yet, check 
			indiv[member].quarantine = 1; // isolate this person
			stateNumber[4]--;
			stateNumber[7]++; // increase a number of isolated people by 1
			
			// clear waiting information of PCR testing of this isolated person
			indiv[member].waitingResult = 0; // relese the person from waiting 
			indiv[member].waitingDays = 0; // clear waiting days (this will be done when they have test, erase the line dose not affect the simulation results)
			}
		}
	} // member
}


void disclosurePCRresult(int stateNumber[], INDIV indiv[]){
	int member, state;
	
	// check if there are individuals waiting for test results
	for(member=0; member< MEMBER; member++){
		if(indiv[member].quarantine == 0){ // if the person has not isolated yet, check 
			if(indiv[member].waitingResult == 1){
				if(indiv[member].waitingDays == REG_READ_TIME){ //see test results
					if(indiv[member].testResult == 1){ // if the person is PCR positive
						// isolate the individual
						indiv[member].quarantine = 1;
						indiv[member].quarantineDays = 0;
						// the person's state
						state = indiv[member].state;
						stateNumber[state]--;
						stateNumber[7]++; // isolation 
					}
					// clear waiting information, for all members
					indiv[member].waitingResult = 0; // relese the person from waiting 
					indiv[member].waitingDays = 0; // clear waiting days (this will be done when they have test, erase the line dose not affect the simulation results)
				} //if(indiv[member].waitingDays == REG_READ_TIME){
			}	
		}
	}
}


int main(){
	int i, j, k, member, week, day, rep, sw, so, dayBegin, whatDay, lastPCR;
	int scenario, state, numInfected, totalNumInfects, dayInfectionCease, totalDayInfectionCease;
	int bp, gameCount, infectedInGame, totalInfectedInGame, quarantineOfTheWeek, massInfection, totalMassInfection;
	int stateNumber[STATES];
	double rnd;
	double beta, gamma, rho, sigma, eta;
	double PCRSTV[STATES], antigenSTV[STATES];
	INDIV indiv[MEMBER];
  
	setRandomSeed();
	
	// parameters
	
	//average duration as E
	sw = 3; //wile type
	so = 1; // omicron
	
	sigma = DELTA * 1.0 / (double)sw; // rate from E(1) to P1, sw for wild type, so for omicron
	rho = DELTA * 1.0 / 1.0; //Ia1 and Ia2 last 1 day
	gamma = DELTA * 1.0 / 7.0; // recovery in 7 days
	eta = 0.54;
	beta = FREQ * DELTA * R_0 / 9.0; 
	
	
	/*testing scenarios
	scenario = 0; daily symptom only
	scenario = 1; bi-weekly PCR
	scenario = 2; weekly PCR
	scenario = 3; twice antigen in a week with 35% relative sensitivity
	scenario = 4; twice antigen in a week with 50% relative sensitivity
	scenario = 5; twice antigen in a week with 75% relative sensitivity
	*/
	for(scenario = 0; scenario<6; scenario++){
		//set sensitivity of the PCR testing
		setPCRSensitivity(PCRSTV);
		setAntigenSensitivity(antigenSTV, PCRSTV, scenario);
	
	
		totalNumInfects = 0;
		totalDayInfectionCease = 0;
		infectedInGame = 0;
		gameCount = 0;
		totalMassInfection = 0;
		
		
		for(rep=0; rep<10000; rep++){
			// counters for measure items
			bp = 0; // reset tag for break
			quarantineOfTheWeek = 0;
			massInfection = 0;
			
			// initialization
			initializePopulation(stateNumber, indiv);
			// a day of week when new E arize
			dayBegin = (int)(7.0 *urand()); //0: Saturday, 1: Sunday,..., 6: Friday
			// When was the last PCR, 0: two weeks ago, 1: last week
			lastPCR = (int)(2.0 * urand());
			// end initialization
		
			// make one E individual
			indiv[0].state = 1;
			stateNumber[0]--;
			stateNumber[1]++;
		
		
			for(week=0; week<38; week++){ // simulation length is 6 weeks
				for(day = dayBegin; day<7+dayBegin; day++){
					//What day is it today?
					whatDay = day%7; //0: Saturday, 1: Sunday,..., 6: Friday
				
					// increment  waiting day for test results
					for(member=0; member<MEMBER; member++){
						if(indiv[member].waitingResult == 1) indiv[member].waitingDays++;
						if(indiv[member].quarantine == 1) indiv[member].quarantineDays++;
					}
			
					// daily symptom check
					dailySymptomCheck(stateNumber, indiv);
					
					// see the results of PCR testing
					disclosurePCRresult(stateNumber, indiv);
					

				// folk by the testing scenario
				switch(scenario) {
				  case 0: // do nothing for testing
					break;
				  case 1: //bi-weekly PCR testing
					if(whatDay == 6){ // if it is Friday today, 
						if(week%2 == lastPCR){ // and if no PCR testing last week
							doTest(indiv, PCRSTV);
						}
					}
					break;
				  case 2: //weekly PCR
					if(whatDay == 6){ // if it is Friday today, 
						doTest(indiv, PCRSTV);
					}
					break;
				  case 3: //twice antigen in a week with 35% relative sensitivity
					if(whatDay == 3 || whatDay == 6){ // if it is Tuesday or Friday
						doAntigenTest(stateNumber, indiv, antigenSTV);
					}
					break;
				  case 4: //twice antigen in a week with 50% relative sensitivity
					if(whatDay == 3 || whatDay == 6){ // if it is Tuesday or Friday
						doAntigenTest(stateNumber, indiv, antigenSTV);
					}
					break;
				  case 5: //twice antigen in a week with 70% relative sensitivity
					if(whatDay == 3 || whatDay == 6){ // if it is Tuesday or Friday
						doAntigenTest(stateNumber, indiv, antigenSTV);
					}
					break;
				} // switch
					
					
					// some for statistics
					if(whatDay == 0){
						//play a game
						gameCount++;
						infectedInGame += (stateNumber[2] + stateNumber[3] + stateNumber[4] + stateNumber[5]);
						//printf("%d %d\n", gameCount, infectedInGame);
					}
					
					
					// proceed infection for one day
					infections_in_a_day(stateNumber, indiv, beta, gamma, rho, sigma, eta);
					
					// 
					if(stateNumber[1] + stateNumber[2] + stateNumber[3] + stateNumber[4] + stateNumber[5] == 0){
						// cease infection
						dayInfectionCease = (day-dayBegin) + 7*week;
						bp=1; // if there are no infected individuals, quit simulation
						break;
					}
				
				} // day
				
				// check for mass infection
				if(massInfection == 0) {
					if(stateNumber[7]-quarantineOfTheWeek > 4){ //mass infection occurs
						massInfection = 1;
					}
					quarantineOfTheWeek = stateNumber[7];
				}
				
				if(bp==1){
					if(stateNumber[7] == 0) zeroFind++;
					break;
				}

		} // week
		
			numInfected = MEMBER-stateNumber[0];
			totalNumInfects += numInfected;
			totalDayInfectionCease += dayInfectionCease;
			totalMassInfection += massInfection;
		
		//printf("infected number = %d\n", numInfected);
		
	} // rep
	
		printf("%g %g %g %g\n", (double)totalNumInfects/10000.0, (double)totalDayInfectionCease/10000.0, (double)infectedInGame/(double)gameCount, (double)totalMassInfection/10000.0);
	}
	
	return 0;
}
