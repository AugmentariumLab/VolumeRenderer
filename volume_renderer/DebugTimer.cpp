// #include "header.h"
#include "DebugTimer.h"

map<string, SingleTimePoint> DebugTimer::Map;

void DebugTimer::Begin(int num_timers, string label) {
#ifdef DEBUG_TIMER_ON
	Map[label].startTime = chrono::steady_clock::now();
	Map[label].numTimers = num_timers;
#endif // DEBUG_TIMER_ON
}
void DebugTimer::End(string label) {
#ifdef DEBUG_TIMER_ON
	SingleTimePoint* t = &Map[label];
	auto end = chrono::steady_clock::now();
	auto diff = end - t->startTime;
	t->totalTime += chrono::duration <double, milli>(diff).count();

	double ms;
	double fps;
	if (t->counter == t->numTimers - 1)
	{
		ms = t->totalTime / t->numTimers;
		fps = 1.0 / (ms / 1000.0);
		cout << "[" << label << "]\t" << t->totalTime / t->numTimers << " ms " << " (" << fps << " fps)" << endl;
		t->counter = 0;
		t->totalTime = 0;
	}
	else
		t->counter++;

	//cout << "[" << label << "]\t" << chrono::duration <double, milli>(diff).count() << " ms" << endl;
	//cout << counter << endl;
#endif // DEBUG_TIMER_ON
}