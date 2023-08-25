#include <gmpxx.h>
#include <assert.h>
#include <iostream>
#include <cstdio>
#include <random>
#include "include/point.h"
#include "include/curve.h"
#include "include/mpz_math.h"
#include "include/peer.h"

#include <arpa/inet.h>
#include <sys/socket.h>
#include <sys/ioctl.h>
#include <net/if.h>
#include <unistd.h>
#include <sodium.h>
#include <cstring>
#include <sstream>
#include <fstream>
#include <x86intrin.h>

#define cpucycles(cycles) cycles = __rdtsc()
#define cpucycles_reset() cpucycles_sum = 0
#define cpucycles_start() cpucycles(cpucycles_before)
#define cpucycles_stop()                                 \
  do                                                     \
  {                                                      \
    cpucycles(cpucycles_after);                          \
    cpucycles_sum += cpucycles_after - cpucycles_before; \
  } while (0)

#define cpucycles_result() cpucycles_sum

unsigned long long cpucycles_before, cpucycles_after, cpucycles_sum;

using namespace std;

#define SERVER_IP "192.168.10.1"    

std::string getMACAddress(const char* interfaceName, int fd) {
    struct ifreq ifr;

    std::strcpy(ifr.ifr_name, interfaceName);

    if (ioctl(fd, SIOCGIFHWADDR, &ifr) == -1) {
        std::cerr << "Error getting MAC address" << std::endl;
        return "";
    }

    char macAddress[18];
    std::snprintf(macAddress, sizeof(macAddress), "%.2X:%.2X:%.2X:%.2X:%.2X:%.2X",
        static_cast<unsigned char>(ifr.ifr_hwaddr.sa_data[0]),
        static_cast<unsigned char>(ifr.ifr_hwaddr.sa_data[1]),
        static_cast<unsigned char>(ifr.ifr_hwaddr.sa_data[2]),
        static_cast<unsigned char>(ifr.ifr_hwaddr.sa_data[3]),
        static_cast<unsigned char>(ifr.ifr_hwaddr.sa_data[4]),
        static_cast<unsigned char>(ifr.ifr_hwaddr.sa_data[5]));
    
    return macAddress;
}

int main(int argc, char* argv[]) {
    srand(time(NULL));

    if (sodium_init() != 0) {
        printf("panic! the library couldn't be initialized; it is not safe to use");
        return 1;
    }

    /**
    mpz_t a, b, c, d;
    mpz_init(a);
    mpz_init(b);
    mpz_init(c);
    mpz_init(d);
    //计算2的1000次方
    mpz_init_set_ui(a, 2);
    mpz_pow_ui(c, a, 1000);
    gmp_printf("c = %Zd\n", c);

    //计算12345678900987654321*98765432100123456789
    mpz_init_set_str(b, "41577311594322086732010735089046630069560194285817217577730709794461777884595", 10); //10进制
    mpz_init_set_str(c, "76884956397045344220809746629001649093037950200943055203735601445031516197751", 10);
    mpz_mul(d, b, c);
    gmp_printf("d = %Zd\n", d);

    //计算除法
    // mpz_init_set_str(c, "9", 10);
    // mpz_init_set_str(b, "2", 10);//10进制
    mpz_fdiv_q(d, c, b); //向下整除
    gmp_printf("mpz_fdiv_q:d = %Zd\n", d);

    legendre(d, b, c);
    gmp_printf("legendre:d = %Zd\n", d);

    tonelli_shanks(d, b, c);
    gmp_printf("tonelli_shanks:d = %Zd\n", d);

    EllipticCurve ec(
        "56698187605326110043627228396178346077120614539475214109386828188763884139993",
        "17577232497321838841075697789794520262950426058923084567046852300633325438902",
        "76884956397045344220809746629001649093037950200943055203735601445031516197751", 10);
    mpz_t x;

    mpz_init_set_str(x, "44636889021951429617878751535492807962563326887208848194489661399766546537096", 10);
    ec.curve_equation(d, x);
    gmp_printf("curve_equation :d = %Zd\n", d);
    mpz_clear(x);

    bool is_QR;
    is_QR = ec.is_quadratic_residue(d);
    gmp_printf("is_quadratic_residue : d = %d\n", is_QR);

    mpz_init_set_str(x, "1128969709662994286855259741431387894916776047019745036539814346559507564901", 10);
    is_QR = ec.is_quadratic_residue(x);
    gmp_printf("is_quadratic_residue :d = %d\n", is_QR);
    mpz_clear(x);

    Point p(
        "11245091364210913538031558554069068295821326089979728553744121109670093869553",
        "58227085528211006908619226448297745234381220210201339909209560974723630792484",
        10);
    bool is_on_curve = ec.is_point_on_curve(p);
    gmp_printf("is_point_on_curve :d = %d\n", is_on_curve);

    mpz_init_set_str(x, "31585362742042248681507141308670880101726446859422047342877633324024259170673", 10);
    ec.inv_mod_p(d, x);
    gmp_printf("inv_mod_p :d = %Zd\n", d);

    Point p1(
        "63463127379121579880278198109336113474044035466585342005071192132059818598209",
        "16558137935435870897706968960488825780425063076841601361802226641039654238867",
        10);
    Point p2;
    ec.ec_point_inv(p2, p1);
    gmp_printf("ec_point_inv :x = %Zd , y = %Zd \n", p2.x, p2.y);

    Point p3(
        "44798834989657919863089174096808959741372081687640931490811108075146543214623",
        " 32552174646444655060134058623257090394879394920791321656975553975238125448510",
        10);
    Point p4;
    ec.ec_point_double_add(p4, p3);
    gmp_printf("ec_point_double_add :x = %Zd , y = %Zd \n", p4.x, p4.y);

    p1.set_x("19255648800734979802927059587974852932026844473492005743471589635415939677617", 10);
    p1.set_y("5347230032267046714807334904413528386238164048132981468634421549716907951517", 10);
    p2.set_x("53966716554608051077465919234931511423879734331111662366752677474299058504440", 10);
    p2.set_y("8117559023792801597289705548012874557664944241595251311237145384718921880283", 10);
    ec.ec_point_add(p4, p1, p2);
    gmp_printf("ec_point_add :x = %Zd , y = %Zd \n", p4.x, p4.y);

    p3.set_x("72639184339425654576555087557111920487359235746629312591395173225124908193763", 10);
    p3.set_y("15224464906049663749133738389570345770876661096229951121126739762545982893453", 10);
    mpz_init_set_str(x, "30365056390190848852943231209253184370325221626750165764764158042136144461840", 10);
    ec.ec_point_scalar_mul(p4, x, p3);
    mpz_clear(x);
    gmp_printf("ec_point_scalar_mul :x = %Zd , y = %Zd \n", p4.x, p4.y);

    p3.set_x("1104832230923804541019855326938339287498267286713819743046729009762879975423", 10);
    p3.set_y("52921860843340567879358341355783097105123683025468323692090409904865326996090", 10);
    mpz_init_set_str(x, "25568832556651077009177958848123513550279603375538798880147678520829252774745", 10);
    ec.ec_point_scalar_mul(p4, x, p3);
    mpz_clear(x);
    gmp_printf("ec_point_scalar_mul :x = %Zd , y = %Zd \n", p4.x, p4.y);
    */

    int server_fd, new_socket, valread, bytes_sent;
    int opt = 1;

    if (argc != 3) {
        printf("Usage: %s <LISTENING PORT> <PEER NAME>\n", argv[0]);
        return 1;
    }

    int PORT = atoi(argv[1]);
    char *SERVER_NAME = argv[2];

    struct sockaddr_in address;
    int addrlen = sizeof(address);
    address.sin_family = AF_INET;
    address.sin_port = htons(PORT);
    address.sin_addr.s_addr = INADDR_ANY;

    printf("[SOCK] Listening for %s machine at port %i\n", SERVER_NAME, PORT);

    // Creating socket file descriptor
    if ((server_fd = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0) {
        std::cerr << "Error opening socket" << std::endl;
        exit(EXIT_FAILURE);
    }
    // Forcefully attaching socket to the port 8080 to SERVER_IP
    if (inet_pton(AF_INET, SERVER_IP, &address.sin_addr) <= 0) {
        std::cerr << "Invalid address/ Address not supported" << std::endl;
        return -1;
    }
    // Forcefully attaching socket
    if (setsockopt(server_fd, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt))) {
        std::cerr << "setsockopt" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (bind(server_fd, (struct sockaddr *)&address, sizeof(address)) < 0) {
        std::cerr << "Bind failed" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (listen(server_fd, 3) < 0) {
        std::cerr << "Error listening" << std::endl;
        exit(EXIT_FAILURE);
    }
    if ((new_socket = accept(server_fd, (struct sockaddr *)&address, (socklen_t *)&addrlen)) < 0) {
        std::cerr << "Error accepting" << std::endl;
        exit(EXIT_FAILURE);
    }

    // std::cout << "=========== Dragonfly Key Exchange ==========" << std::endl;
    std::cout << "[SOCK] Starting up on " << SERVER_IP << ":" << PORT << std::endl;
    
    gmp_printf("---------------------------------------------------\n");

    std::string own_mac = getMACAddress("eth1", server_fd); // Change to your interface name
    Peer ap("abc1238", own_mac, "AP");

    // Start timer
    struct timespec receive_mac, sending_mac, receive_scalar, sending_scalar, receive_token, sending_token;
    struct timespec begin_dragonfly_wall, begin_dragonfly_cpu, end_dragonfly_wall, end_dragonfly_cpu;
    clock_gettime(CLOCK_REALTIME, &begin_dragonfly_wall);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin_dragonfly_cpu);

    uint8_t sta_mac[18];
    valread = recv(new_socket, sta_mac, sizeof(sta_mac), 0);
    clock_gettime(CLOCK_REALTIME, &receive_mac);
    std::cout << "[AP] Received Data Size " << valread << " Bytes." << std::endl;
    std::cout << "[AP] STA MAC Address received." << std::endl;

    bytes_sent = send(new_socket, own_mac.c_str(), 18, 0);
    clock_gettime(CLOCK_REALTIME, &sending_mac);
    std::cout << "\n[AP] Transmitted Data Size " << bytes_sent << " Bytes." << std::endl;
    std::cout << "[AP] AP MAC Address sent." << std::endl;

    double MAC_TIME = (sending_mac.tv_sec - receive_mac.tv_sec) + (sending_mac.tv_nsec - receive_mac.tv_nsec) / 1000000000.0 * 1000.0;

    std::cout << "\n[AP] Own MAC Address: " << own_mac << std::endl;
    std::cout << "[AP] Peer MAC Address: " << sta_mac << std::endl;

    gmp_printf("---------------------------------------------------\n");
    gmp_printf("Starting hunting and pecking to derive PE...\n\n");
    std::string sta_string_mac(reinterpret_cast<char*>(sta_mac), 18);

    cpucycles_reset();
    cpucycles_start();
    ap.initiate(sta_string_mac);
    cpucycles_stop();
    unsigned int sta_initiate_cycles = cpucycles_result();

    gmp_printf("---------------------------------------------------\n");
    gmp_printf("Starting dragonfly commit exchange...\n\n");

    // AP running commit exchange to generate scalar & element
    cpucycles_reset();
    cpucycles_start();
    ap.commit_exchange();
    cpucycles_stop();
    unsigned int sta_commit_cycles = cpucycles_result();

    // Putting Scalar and Element into stream to be sent to peer
    std::ostringstream scalar_element_stream;

    // Stream the scalar and element
    scalar_element_stream << mpz_get_str(nullptr, 10, ap.scalar) << "\n";
    scalar_element_stream << mpz_get_str(nullptr, 10, ap.element.x) << "\n";
    scalar_element_stream << mpz_get_str(nullptr, 10, ap.element.y);
    std::string scalar_complete_stream = scalar_element_stream.str();

    // Receiving peer Scalar & Element
    uint8_t peer_scalarElement[512];
    valread = recv(new_socket, peer_scalarElement, sizeof(peer_scalarElement), 0);
    clock_gettime(CLOCK_REALTIME, &receive_scalar);
    std::cout << "\n[AP] Received Data Size " << valread << " Bytes." << std::endl;
    std::cout << "[AP] STA Scalar Element received." << std::endl;

    // Sending Scalar & Element to peer
    bytes_sent = send(new_socket, scalar_complete_stream.c_str(), 512, 0);
    clock_gettime(CLOCK_REALTIME, &sending_scalar);
    std::cout << "\n[AP] Transmitted Data Size " << bytes_sent << " Bytes." << std::endl;
    std::cout << "[AP] AP Scalar Element sent." << std::endl;

    double SCALAR_TIME = (sending_scalar.tv_sec - receive_scalar.tv_sec) + (sending_scalar.tv_nsec - receive_scalar.tv_nsec) / 1000000000.0 * 1000.0;

    // Splitting received Scalar and Element from stream
    mpz_t scalar_sta, element_x_sta, element_y_sta;
    mpz_init(scalar_sta);
    mpz_init(element_x_sta);
    mpz_init(element_y_sta);

    std::istringstream stream(std::string(reinterpret_cast<char*>(peer_scalarElement)));
    std::string received_scalar, received_elementX, received_elementY;
    while (std::getline(stream, received_scalar) && std::getline(stream, received_elementX) && std::getline(stream, received_elementY)) {
        mpz_set_str(scalar_sta, received_scalar.c_str(), 10);
        mpz_set_str(element_x_sta, received_elementX.c_str(), 10);
        mpz_set_str(element_y_sta, received_elementY.c_str(), 10);
    }

    // Create a Point object using the parsed coordinates
    Point element_sta(element_x_sta, element_y_sta);

    gmp_printf("---------------------------------------------------\n");
    gmp_printf("Computing shared secret...\n\n");

    cpucycles_reset();
    cpucycles_start();
    std::string ap_token = ap.compute_shared_secret(element_sta, scalar_sta, sta_string_mac);
    cpucycles_stop();
    unsigned int sta_ss_cycles = cpucycles_result();

    mpz_clear(element_x_sta);
    mpz_clear(element_y_sta);
    
    gmp_printf("---------------------------------------------------\n");
    gmp_printf("Confirm Exchange...\n");

    uint8_t sta_token[32];
    valread = recv(new_socket, sta_token, 32, 0);
    clock_gettime(CLOCK_REALTIME, &receive_token);
    std::cout << "\n[AP] Received Data Size " << valread << " Bytes." << std::endl;
    std::cout << "[AP] STA token received." << std::endl;

    bytes_sent = send(new_socket, ap_token.c_str(), 32, 0);
    clock_gettime(CLOCK_REALTIME, &sending_token);
    std::cout << "\n[AP] Transmitted Data Size " << bytes_sent << " Bytes." << std::endl;
    std::cout << "[AP] AP token sent." << std::endl;

    double TOKEN_TIME = (sending_token.tv_sec - receive_token.tv_sec) + (sending_token.tv_nsec - receive_token.tv_nsec) / 1000000000.0 * 1000.0;

    std::string sta_token_string(reinterpret_cast<char*>(sta_token), 32);

    cpucycles_reset();
    cpucycles_start();
    std::string pmk = ap.confirm_exchange(sta_token_string);
    cpucycles_stop();
    unsigned int sta_confirm_cycles = cpucycles_result();

    // End timer
    clock_gettime(CLOCK_REALTIME, &end_dragonfly_wall);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_dragonfly_cpu);

    string file_ex = "_pmk.key";
    string pmk_filename = SERVER_NAME + file_ex;

    // write pmk key to file
    ofstream pmkfile(pmk_filename);
    pmkfile << pmk;
    pmkfile.close();

    close(new_socket);

    double dragonfly_wall = (end_dragonfly_wall.tv_sec - begin_dragonfly_wall.tv_sec) + (end_dragonfly_wall.tv_nsec - begin_dragonfly_wall.tv_nsec) / 1000000000.0 * 1000.0;
    double dragonfly_cpu = (end_dragonfly_cpu.tv_sec - begin_dragonfly_cpu.tv_sec) + (end_dragonfly_cpu.tv_nsec - begin_dragonfly_cpu.tv_nsec) / 1000000000.0 * 1000.0;

    std::cout << "\nTime taken for MAC Address: " << MAC_TIME << std::endl;
    std::cout << "Time taken for Commit Exchange: " << SCALAR_TIME << std::endl;
    std::cout << "Time taken for Confirm Exchange: " << TOKEN_TIME << std::endl;

    std::cout << "\nTotal CPU Time: " << dragonfly_cpu << " ms" << std::endl;
    std::cout << "Total WALL Time: " << dragonfly_wall << " ms" << std::endl;

    std::cout << "\nDragonfly Initaite: " << sta_initiate_cycles << " CPU Cycles" << std::endl;
    std::cout << "Dragonfly Commit: " << sta_commit_cycles << " CPU Cycles" << std::endl;
    std::cout << "Compute Shared Secret: " << sta_ss_cycles << " CPU Cycles" << std::endl;
    std::cout << "Dragonfly Confirm: " << sta_confirm_cycles << " CPU Cycles" << std::endl;
    std::cout << "Total: " << sta_initiate_cycles + sta_commit_cycles + sta_ss_cycles + sta_confirm_cycles << " CPU Cycles" << std::endl;

    /**   
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(c);
    mpz_clear(d);
    */
    return 0;
}
