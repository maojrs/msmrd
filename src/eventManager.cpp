//
// Created by maojrs on 4/26/19.
//
#include "eventManager.hpp"

namespace msmrd {
    /**
     *  Implementation of class to manage events (mainly transitions and/or reactions) in MSM/RD algorithm.
     */

    /* Adds an event to the event list. This requires specifying the wait time until event happens, the end state
     * for the event/transition, the indexes (in partList) of the particles involved (part1Index < part2Index) and
     * a string: "in" or "out" to show it is transitioning into a bound state or out of it. */
    void eventManager::addEvent(double waitTime, int endState, int part1Index, int part2Index, std::string inORout) {
        std::array<int, 2> partIndexes = {part1Index, part2Index};
        auto event = std::make_tuple(waitTime, endState, partIndexes, inORout);
        eventList.push_back(event);
    }


    /*
     * Routines to remove events, note the different possibilities: one using index in the event list, the second one
     * using the event itself and the third one using a list of events to be erased at once.
     */
    void eventManager::removeEvent(int index){
        eventList.erase(eventList.begin() + index);
    }

    void eventManager::removeEvent(std::tuple<double, int, std::array<int,2>, std::string> event){
        auto currentSize = eventList.size();
        std::remove(eventList.begin(), eventList.end(), event);
        eventList.resize(currentSize - 1);
    }

    void eventManager::removeEvents(std::vector<std::tuple<double, int, std::array<int,2>, std::string>> events) {
        // Not clear yet hot to do this, perhaps: std::remove(eventList.begin(), eventList.end(), events);
    }


    // Returns event in eventList corresponding to the index provided
    std::tuple<double, int, std::array<int,2>, std::string> eventManager::getEvent(int index) {
        return eventList[index];
    }


    // Returns time of event in eventList corresponding to the index provided
    double eventManager::getEventTime(int index) {
        auto event = getEvent(index);
        return std::get<0>(event);
    }


    void eventManager::advanceTime(double timeStep) {
        double newTime;
        for (auto &event : eventList) {
            newTime = std::get<0>(event) - timeStep;
            std::get<0>(event) = newTime;
        }
    }


    // Sorts event list in ascending order using the value of the first entry in the tuple (time).
    void eventManager::sortAscending() {
        std::sort(eventList.begin(), eventList.end());
    }

    // Sorts event list in descending order using the value of the first entry in the tuple (time).
    void eventManager::sortDescending() {
        // Lambda function to define sorting
        auto sortdesc = [](const std::tuple<double, int, std::array<int,2>, std::string>& a,
                           const std::tuple<double, int, std::array<int,2>, std::string>& b) {
            return (std::get<0>(a) > std::get<0>(b));
        };
        std::sort(eventList.begin(), eventList.end(), sortdesc);
    }
}
