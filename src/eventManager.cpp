//
// Created by maojrs on 4/26/19.
//
#include "eventManager.hpp"

namespace msmrd {
    /**
     *  Implementation of class to manage events (mainly transitions and/or reactions) in MSM/RD algorithm.
     */

    /* Adds an event to the event dictionary. This requires specifying the wait time until event happens, the end state
     * for the event/transition, the indexes (in partList) of the particles involved (part1Index < part2Index) and
     * a string: "in" or "out". */
    void eventManager::addEvent(double waitTime, int endState, int part1Index, int part2Index,
                                int originState, std::string inORout) {
        if ((inORout != "in") and (inORout != "out") and (inORout != "inside")) {
            std::range_error("Events can only take 'in', 'out' or 'inside' strings as last argument");
        }
        // Create key and new event
        std::string eventKey = std::to_string(part1Index) + "--" + std::to_string(part2Index);
        event thisEvent = {.waitTime = waitTime, .endState = endState, .part1Index = part1Index,
                           .part2Index = part2Index, originState = originState, .inORout = inORout};
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

    // Advances time of all events
    void eventManager::advanceTime(double timeStep) {
        double newTime;
        for (auto &thisEvent : eventDictionary) {
            thisEvent.second.waitTime = thisEvent.second.waitTime - timeStep;
        }
    }

    // Outputs event log into file
    void eventManager::outputEventLog(int timeIteration, std::string baseFilename) {
        std::ofstream outputfile;
        outputfile.open (baseFilename + "_t" + std::to_string(timeIteration)+ ".txt");
        for (auto &thisEvent : eventDictionary) {
            auto value = thisEvent.second;
            outputfile << value.waitTime << " " << value.endState << " " << value.part1Index <<
                       " " << value.part2Index << " " << value.originState << " " << value.inORout << " \n";
        }
        outputfile.close();
    }

}
