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
     * a string eventType: "binding", "unbinding", "bound2bound", "transition2transition" and "empty". */
    void eventManager::addEvent(double waitTime, int part1Index, int part2Index,
                                int originState, int endState, std::string eventType) {
        if ((eventType != "binding") and (eventType != "unbinding") and
            (eventType != "bound2bound") and (eventType != "transition2transition") and (eventType != "empty")) {
            throw std::invalid_argument("Events can only take 'binding', 'unbinding', 'bound2bound', "
                                        "'transition2transition' or 'empty' strings as last argument");
        }
        // Create key and new event
        std::string eventKey = std::to_string(part1Index) + "--" + std::to_string(part2Index);
        event thisEvent = {.waitTime = waitTime, .part1Index = part1Index, .part2Index = part2Index,
                           .originState = originState, .endState = endState, .eventType = eventType};
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

    /* Removes event associated with particles with indexes: part1Index and part2Index. It is assumed
     * there is only one event associated to a pair of particles and that part1Index < part2Index */
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

    void eventManager::setEventType(std::string newEventType, int part1Index, int part2Index) {
        std::string eventKey = std::to_string(part1Index) + "--" + std::to_string(part2Index);
        auto numKeys = eventDictionary.count(eventKey);
        if (numKeys == 1) {
            eventDictionary[eventKey].eventType = newEventType;
        }
    }

    // Advances time of all events
    void eventManager::advanceTime(double timeStep) {
        double newTime;
        for (auto &thisEvent : eventDictionary) {
            thisEvent.second.waitTime = thisEvent.second.waitTime - timeStep;
        }
    }


    // Writes current events to event logfile
    void eventManager::write2EventLog(int timeIteration) {
        eventLog.push_back(std::to_string(timeIteration));
        if (eventDictionary.size() > 0) {
            for (auto &thisEvent : eventDictionary) {
                auto value = thisEvent.second;
                auto event = std::to_string(value.waitTime) + " " + std::to_string(value.part1Index) +
                         " " + std::to_string(value.part2Index) + " " + std::to_string(value.originState) +
                         " " + std::to_string(value.endState) + " " + value.eventType + " \n";
                eventLog.push_back(event);
            }
        } else{
            eventLog.push_back("No events \n");
        }
    }

    // Prints log into file
    void eventManager::printEventLog(std::string baseFilename) {
            std::ofstream outputfile;
            outputfile.open (baseFilename + ".dat");
            for (auto &thisLine : eventLog) {
                outputfile << thisLine << " \n";
            }
            outputfile.close();
            eventLog.clear();
    }


}
