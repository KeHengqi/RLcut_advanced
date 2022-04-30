from fileinput import close
import os
import re
Strategy = ["in-going", "degree", "DC_speed", "DC_price", "mirror"]
# A complite graph file should be: ./dataset + Graph_file
# this is for tencent server
# Graph_file = ["debug/debug.txt", "facebook/facebook_combined.txt", "gnutella/p2p-Gnutella09.txt", "googleweb/web-Google.txt", "Wiki-Vote/Wiki-Vote.txt"]

# this is for school server
# on school server, the graph file should be at /home/local_graph
Graph_file = ["debug.txt", "it", "orkut", "googleweb/web-Google.txt", "twitter", "livejournal", "soc-LiveJournal1.txt", "uk"]

# this is for tencent server
# Graph_name = ["debug", "facebook", "gnutella", "google", "wiki"]

# this is for school server
Graph_name = ["debug", "it", "orkut", "google", "twitter", "livejournal", "soc-LiveJournal1", "uk"]

# A complite network file should be: ./network + Network_file
# this is for tencent server
# Network_file = ["debug.txt", "amazon", "azure",  "Network_cost.txt", "Network_high.txt", "Network_homo.txt", "Network_medium.txt", "Network.txt"]

# this is for school server
Network_file = ["debug.txt", "amazon", "amazon_high", "amazon_low", "amazon_medium", "azure",  "Network_cost.txt", "Network_high.txt", "Network_low.txt", "Network_medium.txt", "Network.txt"]
# print(DC_number)


# default_sentence = "./main pagerank ./dataset/Wiki-Vote/Wiki-Vote.txt ./network/amazon 8 0.4 100 48 output/wiki 1 train/wiki  1 8 1000000000 1 0 "
default_sentence = "./main "
TUI_OR_SCRIPT = "script "
Application = "pagerank "
# Can add a choice for application, not now

default_sentence += TUI_OR_SCRIPT + Application

print("Welcome to the graph seperation execution script")
print("First choose the graph you want to test: ")

for counter in range(1, len(Graph_name)):
    print(str(counter) + ". " + Graph_name[counter])

# print("1. facebook\n2. gnutella\n3. google\n4. wiki")
while True:
    graph_choice = int(input("Now choose one: "))
    if graph_choice < 0 or graph_choice >= len(Graph_name):
        print("Illegal input! Please input again!")
        print("----------------------------------")
    else:
        break
# for tencent server
# graph_data = "./dataset/" + Graph_file[graph_choice]

# for school server
graph_data = "/home/local_graph/" + Graph_file[graph_choice]
default_sentence += graph_data + " "

print("Now choose the network you want to simulate")
for counter in range(1, len(Network_file)):
    print(str(counter) + ". " + Network_file[counter])

# print("1. amazon\n2. azure\n3. Network_cost\n4. Network_high\n5. Network_homo\n6. Network_medium\n7. Network")
while True:
    network_choice = int(input("Now choose one: "))
    if network_choice < 0 or network_choice > 7:
        print("Illegal input! Please input again!")
        print("----------------------------------")
    else:
        break

# get the DC number
chosen_network_file = "./network/" + Network_file[int(network_choice)]
network = open(chosen_network_file, "r")
DC_number = 0
while network.readline():
    DC_number += 1

default_sentence += chosen_network_file + " " + str(DC_number) + " "

budget = 0.4
theta = 100
thread_num = 48
randpro = 1
bsp = 1
batchsize = 8
left_time = 1000000000
n = 1
ginger_cost = 0
default_sentence += str(budget) + " " + str(theta) + " " + str(thread_num) + " "


print("Here are the option you can choose")
print("1. in-going\n2. degree\n3. DC_speed\n4. DC_price\n5. mirror")
option = input("Now choose one to start the script: ")
character_para = Strategy[int(option) - 1]
# if int(option) == 1:
#     character_para += " n"
#     print(character_para)

output_folder = "./output/" + Strategy[int(option) - 1] + "/" + Graph_name[graph_choice] + "/"
if os.path.isdir(output_folder) == False:
    path = os.path.join(output_folder)
    os.makedirs(path)

trained_folder = "./train/" + Strategy[int(option) - 1] + "/" + Graph_name[graph_choice] + "/"
if os.path.isdir(trained_folder) == False:
    path = os.path.join(trained_folder)
    os.makedirs(path)
test_times = 5

all_output_file = list()
all_trained_file = list()

for counter in range(1, 11):
    ratio = counter / 10
    all_output_file.clear()
    all_trained_file.clear()
    output_ratio_folder = output_folder + str(ratio) + "/"
    if os.path.isdir(output_ratio_folder) == False:
        path = os.path.join(output_ratio_folder)
        os.makedirs(path)
    trained_ratio_folder = trained_folder + str(ratio) + "/"
    if os.path.isdir(trained_ratio_folder) == False:
        path = os.path.join(trained_ratio_folder)
        os.makedirs(path)
    # print(ratio)
    output_file_common = output_ratio_folder + Graph_name[graph_choice] + "-"
    trained_file_common = trained_ratio_folder + Graph_name[graph_choice] + "-"
    for file_no in range(0, test_times):
        output_file = output_file_common + str(ratio) + "-" + str(file_no + 1) + ".csv"
        all_output_file.append(output_file)
        trained_file = trained_file_common + str(ratio) + "-" + str(file_no + 1) + ".txt"
        all_trained_file.append(trained_file)

        character = character_para + " " + str(ratio) + " "
        if int(option) == 1:
            character += "n"
        # print(character_para)
        running_sentence = default_sentence + output_file + " " + str(randpro) + " " + trained_file + " " + str(bsp) + " " + str(batchsize) + " " + str(left_time) + " " + str(n) + " " + str(ginger_cost) + " "
        running_sentence += character
        os.system(running_sentence)

    output_file_average = output_ratio_folder + Graph_name[graph_choice] + "-" + str(ratio) + "-average.csv"
    time_cost_list = list()
    price_cost_list = list()
    overhead_list = list()
    line_number_list = list()
    line_number_list.append(0)

    for distributed_file in all_output_file:
        dfopen = open(distributed_file, "r")
        line_number = 0
        line_content = dfopen.readline()
        while line_content:
            #print(line_content)
            line_content = dfopen.readline()
            line_number += 1
            # print(line_content)
            line_content_token_list = re.split(r'[, \n]', line_content)
            if len(line_content_token_list) == 1:
                break
            # print(line_content_token_list)
            # print(line_content_token_list[1])
            # print(line_content_token_list[2])
            # print(line_content_token_list[3])
            if line_content_token_list[1] != '':
                time_cost_list.append(float(line_content_token_list[1]))
            else:
                time_cost_list.append(0)

            if line_content_token_list[2] != '':
                price_cost_list.append(float(line_content_token_list[2]))
            else:
                price_cost_list.append(0)
            # time_cost_list.append(float(line_content_token_list[1]))
            # price_cost_list.append(float(line_content_token_list[2]))
            if line_content_token_list[3] != '':
                overhead_list.append(float(line_content_token_list[3]))
            else:
                overhead_list.append(0)
        line_number -= 1
        line_number_list.append(line_number)
    
    print(line_number_list)
    print(len(line_number_list))
    out_file_avg_open = open(output_file_average, "w")
    prefix = Graph_name[graph_choice] + "-"+ str(ratio) + "-"
    out_file_avg_open.write(prefix + "iteration," + prefix + "time_cost," + prefix + "price_cost," + prefix + "overhead\n")
    # Attention: Suppose all files have same iteration
    for distributed_number in range(0, line_number_list[1]):
        time_cost_avg = 0
        price_cost_avg = 0
        overhead_avg = 0
        index = distributed_number
        for file_number in range(0, len(line_number_list) - 1):
            time_cost_avg += time_cost_list[index + line_number_list[file_number]]
            price_cost_avg += price_cost_list[index + line_number_list[file_number]]
            overhead_avg += overhead_list[index + line_number_list[file_number]]
            index += line_number_list[file_number]
        time_cost_avg /= test_times
        price_cost_avg /= test_times
        overhead_avg /= test_times
        line_data = str(distributed_number) + "," + str(time_cost_avg) + "," + str(price_cost_avg) + ","
        if distributed_number == 0:
            line_data += "\n"
        else:
            line_data += str(overhead_avg) + "\n"
        out_file_avg_open.write(line_data)

    # close(out_file_avg_open)
    #print(running_sentence)
    # run the running sentence

    # if os.system(running_sentence):
    #     print("The script is running well!")
    # else:
    #     print("The script seems to have some problem, please check")

    #print(running_sentence)

""" if os.system():
    print("The script is running well!")
else:
    print("The script seems to have some problem, please check") """