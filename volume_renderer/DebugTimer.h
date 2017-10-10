#pragma once
#include <map>
#include <string>
#include <chrono>
#include <vector>
#include <iostream>
using namespace std;
#define DEBUG_TIMER_ON
//#define TOTAL_TIMERS 250


struct SingleTimePoint {
	chrono::steady_clock::time_point startTime;
	int counter;
	double totalTime;
	int numTimers;
};

class DebugTimer {
private:
	static map<string, SingleTimePoint> Map;
public:
	static void Begin(int num_timers, string label = "-");
	static void End(string label = "-");
};