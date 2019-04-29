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
    public:
        /*
         * @param eventList, list of events/reaction/transisition. Each entry of the list has 4 components, time until
         * event happens, the state to which it will transision (endstate), an array containing the two indexes (in the
         * particleList) of the particles involved in the event and a string, "in" or "out" to show it is transitioning
         * to into a bound state or out of it. The list should be always sorted after being updated by ascending time
         * (first entry).
         */
        std::vector<std::tuple<double, int, std::array<int,2>, std::string>> eventList  = {};

        eventManager() {};

        void addEvent(double waitTime, int endState, int part1Index, int part2Index, std::string inORout);

        void removeEvent(int index);

        void removeEvent(std::tuple<double, int, std::array<int,2>, std::string> event);

        void removeEvents(std::vector<std::tuple<double, int, std::array<int,2>, std::string>> events);

        std::tuple<double, int, std::array<int,2>, std::string> getEvent(int index);

        double getEventTime(int index);

        void advanceTime(double timeStep);

        int getNumEvents() { return eventList.size(); }

        void sort();
    };
}
