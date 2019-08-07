//
// Created by maojrs on 4/26/19.
//

#pragma once
#include <bits/stdc++.h>
#include <tuple>
#include <vector>


namespace msmrd {
    /**
     * Class to manage events (mainly transitions and/or reactions) in MSM/RD algorithm.
     */
    class eventManager {
    private:
        struct event {
            double waitTime;
            int part1Index;
            int part2Index;
            int originState;
            int endState;
            std::string inORout;
        };
        struct event emptyEvent = {.waitTime = std::numeric_limits<double>::infinity(),
               .part1Index = -1, .part2Index = -1, .originState = -1,  .endState = -1, .inORout = "empty"};
    public:
        std::map<std::string, event> eventDictionary;

        /*
         * @struct event, each event consists of 5 variables, time until event happens, the state to which it will
         * transision (endstate), indexes (in the particleList) of the particles involved in the event
         * (part1Index < part2Index always) and a string, "in", "out" or "inside" to show it is transitioning into
         * a bound state, out of it or within bound states, repsectively.
         * @param emptyEvent defines default emptyEvent, returns event with waitTime = infinity.
         * @param eventDictionary, list of events/reaction/transitions, associated with a key corresponding to the
         * two indexes of the particles involved, e.g. "part1Index--part2Index", or more specifically "12--2".
         */

        eventManager() {};

        void addEvent(double waitTime, int part1Index, int part2Index,
                      int originState, int endState, std::string inORout);

        void removeEvent(int part1Index, int part2Index);

        double getEventTime(int part1Index, int part2Index);

        struct event getEvent(int part1Index, int part2Index);

        void advanceTime(double timeStep);

        int getNumEvents() { return eventDictionary.size(); }

        void outputEventLog(int timeIteration, std::string baseFilename);

    };
}
