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
//    void eventManager::addEvent(double waitTime, int endState, int part1Index, int part2Index, std::string inORout) {
//        if ((inORout != "in") and (inORout != "out")) {
//            std::range_error("Events can only take 'in' or 'out' strings as last argument");
//        }
//        std::array<int, 2> partIndexes = {part1Index, part2Index};
//        auto event = std::make_tuple(waitTime, endState, partIndexes, inORout);
//        eventList.push_back(event);
//    }

    /* Adds an event to the event dictionary. This requires specifying the wait time until event happens, the end state
     * for the event/transition, the indexes (in partList) of the particles involved (part1Index < part2Index) and
     * a string: "in" or "out". */
    void eventManager::addEvent(double waitTime, int endState, int part1Index, int part2Index, std::string inORout) {
        if ((inORout != "in") and (inORout != "out")) {
            std::range_error("Events can only take 'in' or 'out' strings as last argument");
        }
        // Create key and new event
        std::string eventKey = std::to_string(part1Index) + "--" + std::to_string(part2Index);
        event thisEvent = {.waitTime = waitTime, .endState = endState, .part1Index = part1Index,
                           .part2Index = part2Index, .inORout = inORout};
        // See if value already exists and if so keep the one with the smallest time. Otherwise insert new event.
        auto search = eventDictionary.find(eventKey);
        // If key and value already assigned in dictionary.
        if (search != eventDictionary.end()) {
            if (search->second.waitTime > waitTime) {
                // Rewrite dictionary value with a given key
                eventDictionary.at(eventKey) = thisEvent;
            }
        } else {
            // Otherwise, insert event into dictionary.
            eventDictionary.insert ( std::pair<std::string, event>(eventKey, thisEvent) );
        }
    }

    void eventManager::removeEvent(int part1Index, int part2Index) {
        std::string eventKey = std::to_string(part1Index) + "--" + std::to_string(part2Index);
        eventDictionary.erase(eventKey);
    }

//    /*
//     * Routines to remove events, note the different possibilities: one using index in the event list, the second one
//     * using the event itself and the third one using a list of events to be erased at once.
//     */
//    void eventManager::removeEvent(int index){
//        eventList.erase(eventList.begin() + index);
//    }

//    void eventManager::removeEvent(std::tuple<double, int, std::array<int,2>, std::string> event){
//        auto currentSize = eventList.size();
//        std::remove(eventList.begin(), eventList.end(), event);
//        eventList.resize(currentSize - 1);
//    }


//    // Returns time of event in eventList corresponding to the index provided
//    double eventManager::getEventTime(int index) {
//        auto event = getEvent(index);
//        return std::get<0>(event);
//    }

    // Returns time of event in eventList corresponding to the index provided
    double eventManager::getEventTime(int part1Index, int part2Index) {
        std::string eventKey = std::to_string(part1Index) + "--" + std::to_string(part2Index);
        auto search = eventDictionary.find(eventKey);
        if (search != eventDictionary.end()) {
            return search->second.waitTime;
        } else {
            return std::numeric_limits<double>::infinity();
        }
    }

    // Returns event in eventList corresponding to the index provided
    struct eventManager::event eventManager::getEvent(int part1Index, int part2Index) {
        std::string eventKey = std::to_string(part1Index) + "--" + std::to_string(part2Index);
        auto search = eventDictionary.find(eventKey);
        if (search != eventDictionary.end()) {
            return search->second;
        } else {
            return emptyEvent;
        }
    }

//    // Returns event in eventList corresponding to the index provided
//    std::tuple<double, int, std::array<int,2>, std::string> eventManager::getEvent(int index) {
//        if (index < eventList.size()) {
//            return eventList[index];
//        } else {
//            return {};
//        }
//    }


    void eventManager::advanceTime(double timeStep) {
        double newTime;
        for (auto &thisEvent : eventDictionary) {
            thisEvent.second.waitTime = thisEvent.second.waitTime - timeStep;
        }
    }

//    void eventManager::advanceTime(double timeStep) {
//        double newTime;
//        for (auto &event : eventList) {
//            newTime = std::get<0>(event) - timeStep;
//            std::get<0>(event) = newTime;
//        }
//    }


//    // Sorts event list in ascending order using the value of the first entry in the tuple (time).
//    void eventManager::sortAscending() {
//        std::sort(eventList.begin(), eventList.end());
//    }
//
//    // Sorts event list in descending order using the value of the first entry in the tuple (time).
//    void eventManager::sortDescending() {
//        // Lambda function to define sorting
//        auto sortdesc = [](const std::tuple<double, int, std::array<int,2>, std::string>& a,
//                           const std::tuple<double, int, std::array<int,2>, std::string>& b) {
//            return (std::get<0>(a) > std::get<0>(b));
//        };
//        std::sort(eventList.begin(), eventList.end(), sortdesc);
//    }
}
